/*------------------------------------------------------------------------------
* postpos.c : post-processing positioning
*
*          Copyright (C) 2007-2016 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/05/08  1.0  new
*           2008/06/16  1.1  support binary inputs
*           2009/01/02  1.2  support new rtk positioing api
*           2009/09/03  1.3  fix bug on combined mode of moving-baseline
*           2009/12/04  1.4  fix bug on obs data buffer overflow
*           2010/07/26  1.5  support ppp-kinematic and ppp-static
*                            support multiple sessions
*                            support sbas positioning
*                            changed api:
*                                postpos()
*                            deleted api:
*                                postposopt()
*           2010/08/16  1.6  fix bug sbas message synchronization (2.4.0_p4)
*           2010/12/09  1.7  support qzss lex and ssr corrections
*           2011/02/07  1.8  fix bug on sbas navigation data conflict
*           2011/03/22  1.9  add function reading g_tec file
*           2011/08/20  1.10 fix bug on freez if solstatic=single and combined
*           2011/09/15  1.11 add function reading stec file
*           2012/02/01  1.12 support keyword expansion of rtcm ssr corrections
*           2013/03/11  1.13 add function reading otl and erp data
*           2014/06/29  1.14 fix problem on overflow of # of satellites
*           2015/03/23  1.15 fix bug on ant type replacement by rinex header
*                            fix bug on combined filter for moving-base mode
*           2015/04/29  1.16 fix bug on reading rtcm ssr corrections
*                            add function to read satellite fcb
*                            add function to read stec and troposphere file
*                            add keyword replacement in dcb, erp and ionos file
*           2015/11/13  1.17 add support of L5 antenna phase center parameters
*                            add *.stec and *.trp file for ppp correction
*           2015/11/26  1.18 support opt->freqopt(disable L2)
*           2016/01/12  1.19 add carrier-phase bias correction by ssr
*           2016/07/31  1.20 fix error message problem in rnx2rtkp
*           2016/08/29  1.21 suppress warnings
*           2016/10/10  1.22 fix bug on identification of file fopt->blq
*           2017/06/13  1.23 add smoother of velocity solution
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define MIN(x,y)    ((x)<(y)?(x):(y))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

#define MAXPRCDAYS  100          /* max days of continuous processing */
#define MAXINFILE   1000         /* max number of input files */

/* constants/global variables ------------------------------------------------*/

static pcvs_t pcvss = { 0 };        /* receiver antenna parameters */
static pcvs_t pcvsr = { 0 };        /* satellite antenna parameters */
static obs_t obss = { 0 };          /* observation data */
static nav_t navs = { 0 };          /* navigation data */
static sbs_t sbss = { 0 };          /* sbas messages */
static lex_t lexs = { 0 };          /* lex messages */
static sta_t stas[MAXRCV];      /* station information */
static int nepoch = 0;            /* number of observation epochs */
static int nitm = 0;            /* number of invalid time marks */
static int iobsu = 0;            /* current rover observation data index */
static int iobsr = 0;            /* current reference observation data index */
static int isbs = 0;            /* current sbas message index */
static int ilex = 0;            /* current lex message index */
static int iitm = 0;            /* current invalid time mark index */
static int revs = 0;            /* analysis direction (0:forward,1:backward) */
static int aborts = 0;            /* abort status */
static sol_t* solf;             /* forward solutions */
static sol_t* solb;             /* backward solutions */
static double* rbf;             /* forward base positions */
static double* rbb;             /* backward base positions */
static int isolf = 0;             /* current forward solutions index */
static int isolb = 0;             /* current backward solutions index */
static char proc_rov[64] = "";   /* rover for current processing */
static char proc_base[64] = "";   /* base station for current processing */
static char rtcm_file[1024] = ""; /* rtcm data file */
static char rtcm_path[1024] = ""; /* rtcm data path */
static gtime_t invalidtm[100] = { {0} };/* invalid time marks */
static rtcm_t rtcm;             /* rtcm control struct */
static FILE* fp_rtcm = NULL;      /* rtcm data file pointer */

/* QHY Defined variables & functions -----------------------------------------*/
#define SATCMP 106

static double lam_1 = CLIGHT / FREQ1_CMP;
static double lam_2 = CLIGHT / FREQ2_CMP;
static double lam_wl = CLIGHT / (FREQ1_CMP - FREQ2_CMP);
static double a_wl = FREQ1_CMP / (FREQ1_CMP - FREQ2_CMP);
static double b_wl = FREQ2_CMP / (FREQ1_CMP - FREQ2_CMP);
static double a_nl = FREQ1_CMP / (FREQ1_CMP + FREQ2_CMP);
static double b_nl = FREQ2_CMP / (FREQ1_CMP + FREQ2_CMP);

static int eph_num;                 /* obss中历元个数 */
static int obs_idx = 0;             /* 记录obs遍历的位置 */
static stec_q stec_list = { 0 };    /* 储存stec */

static int caleph(obs_t* obs);
static double cal_MW(double L1, double L2, double P1, double P2);
static double cal_GF(double L1, double L2, double P1, double P2);
static void sele_satsys(obs_t* obs, prcopt_t* popt, stec_q* stecl);
static void detslp(stec_q* stecl);
static void repslp(double* N_mw, double* N_gf, double* mean_Nmw_p, double* mean_Ngf_p, double* sqr);
static void calSTEC(stec_q* stecl);

/* show message and check break ----------------------------------------------*/
static int checkbrk(const char* format, ...)
{
	va_list arg;
	char buff[1024], * p = buff;
	if (!*format) return showmsg("");

	va_start(arg, format);
	p += vsprintf(p, format, arg);
	va_end(arg);

	if (*proc_rov && *proc_base) sprintf(p, " (%s-%s)", proc_rov, proc_base);
	else if (*proc_rov) sprintf(p, " (%s)", proc_rov);
	else if (*proc_base) sprintf(p, " (%s)", proc_base);

	return showmsg(buff);
}
/* output reference position -------------------------------------------------*/
static void outrpos(FILE* fp, const double* r, const solopt_t* opt)
{
	double pos[3], dms1[3], dms2[3];
	const char* sep = opt->sep;

	trace(3, "outrpos :\n");

	if (opt->posf == SOLF_LLH || opt->posf == SOLF_ENU) {
		ecef2pos(r, pos);
		if (opt->degf) {
			deg2dms(pos[0] * R2D, dms1, 5);
			deg2dms(pos[1] * R2D, dms2, 5);
			fprintf(fp, "%3.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f%s%10.4f",
				dms1[0], sep, dms1[1], sep, dms1[2], sep, dms2[0], sep, dms2[1],
				sep, dms2[2], sep, pos[2]);
		}
		else {
			fprintf(fp, "%13.9f%s%14.9f%s%10.4f", pos[0] * R2D, sep, pos[1] * R2D,
				sep, pos[2]);
		}
	}
	else if (opt->posf == SOLF_XYZ) {
		fprintf(fp, "%14.4f%s%14.4f%s%14.4f", r[0], sep, r[1], sep, r[2]);
	}
}
/* output header -------------------------------------------------------------*/
static void outheader(FILE* fp, char** file, int n, const prcopt_t* popt,
	const solopt_t* sopt)
{
	const char* s1[] = { "GPST","UTC","JST" };
	gtime_t ts, te;
	double t1, t2;
	int i, j, w1, w2;
	char s2[32], s3[32];

	trace(3, "outheader: n=%d\n", n);

	if (sopt->posf == SOLF_NMEA || sopt->posf == SOLF_STAT) {
		return;
	}
	if (sopt->outhead) {
		if (!*sopt->prog) {
			fprintf(fp, "%s program   : RTKLIB ver.%s %s\n", COMMENTH, VER_RTKLIB, PATCH_LEVEL);
		}
		else {
			fprintf(fp, "%s program   : %s\n", COMMENTH, sopt->prog);
		}
		for (i = 0; i < n; i++) {
			fprintf(fp, "%s inp file  : %s\n", COMMENTH, file[i]);
		}
		for (i = 0; i < obss.n; i++)    if (obss.data[i].rcv == 1) break;
		for (j = obss.n - 1; j >= 0; j--) if (obss.data[j].rcv == 1) break;
		if (j < i) { fprintf(fp, "\n%s no rover obs data\n", COMMENTH); return; }
		ts = obss.data[i].time;
		te = obss.data[j].time;
		t1 = time2gpst(ts, &w1);
		t2 = time2gpst(te, &w2);
		if (sopt->times >= 1) ts = gpst2utc(ts);
		if (sopt->times >= 1) te = gpst2utc(te);
		if (sopt->times == 2) ts = timeadd(ts, 9 * 3600.0);
		if (sopt->times == 2) te = timeadd(te, 9 * 3600.0);
		time2str(ts, s2, 1);
		time2str(te, s3, 1);
		fprintf(fp, "%s obs start : %s %s (week%04d %8.1fs)\n", COMMENTH, s2, s1[sopt->times], w1, t1);
		fprintf(fp, "%s obs end   : %s %s (week%04d %8.1fs)\n", COMMENTH, s3, s1[sopt->times], w2, t2);
	}
	if (sopt->outopt) {
		outprcopt(fp, popt);
	}
	if (PMODE_DGPS <= popt->mode && popt->mode <= PMODE_FIXED && popt->mode != PMODE_MOVEB) {
		fprintf(fp, "%s ref pos   :", COMMENTH);
		outrpos(fp, popt->rb, sopt);
		fprintf(fp, "\n");
	}
	if (sopt->outhead || sopt->outopt) fprintf(fp, "%s\n", COMMENTH);

	outsolhead(fp, sopt);
}
/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t* obs, int* i, int rcv)
{
	double tt;
	int n;

	for (; *i < obs->n; (*i)++) if (obs->data[*i].rcv == rcv) break; /* 通过遍历检索接收机号rcv对应的第一个data的位置i */
	for (n = 0; *i + n < obs->n; n++) {
		tt = timediff(obs->data[*i + n].time, obs->data[*i].time);
		if (obs->data[*i + n].rcv != rcv || tt > DTTOL) break; /* 遍历得到当前历元的观测值个数n */
	}
	return n;
}
static int nextobsb(const obs_t* obs, int* i, int rcv)
{
	double tt;
	int n;

	for (; *i >= 0; (*i)--) if (obs->data[*i].rcv == rcv) break;
	for (n = 0; *i - n >= 0; n++) {
		tt = timediff(obs->data[*i - n].time, obs->data[*i].time);
		if (obs->data[*i - n].rcv != rcv || tt < -DTTOL) break;
	}
	return n;
}
/* update rtcm ssr correction ------------------------------------------------*/
static void update_rtcm_ssr(gtime_t time)
{
	char path[1024];
	int i;

	/* open or swap rtcm file */
	reppath(rtcm_file, path, time, "", "");

	if (strcmp(path, rtcm_path)) {
		strcpy(rtcm_path, path);

		if (fp_rtcm) fclose(fp_rtcm);
		fp_rtcm = fopen(path, "rb");
		if (fp_rtcm) {
			rtcm.time = time;
			input_rtcm3f(&rtcm, fp_rtcm);
			trace(2, "rtcm file open: %s\n", path);
		}
	}
	if (!fp_rtcm) return;

	/* read rtcm file until current time */
	while (timediff(rtcm.time, time) < 1E-3) {
		if (input_rtcm3f(&rtcm, fp_rtcm) < -1) break;

		/* update ssr corrections */
		for (i = 0; i < MAXSAT; i++) {
			if (!rtcm.ssr[i].update ||
				rtcm.ssr[i].iod[0] != rtcm.ssr[i].iod[1] ||
				timediff(time, rtcm.ssr[i].t0[0]) < -1E-3) continue;
			navs.ssr[i] = rtcm.ssr[i];
			rtcm.ssr[i].update = 0;
		}
	}
}
/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t* obs, int solq, const prcopt_t* popt)
{
	gtime_t time = { 0 };
	int i, nu, nr, n = 0;

	trace(3, "\ninfunc  : revs=%d iobsu=%d iobsr=%d isbs=%d\n", revs, iobsu, iobsr, isbs);

	if (0 <= iobsu && iobsu < obss.n) {
		settime((time = obss.data[iobsu].time)); /* 设置历元时间 */
		if (checkbrk("processing : %s Q=%d", time_str(time, 0), solq)) {
			aborts = 1; showmsg("aborted"); return -1;
		}
	}
	if (!revs) { /* input forward data */
		if ((nu = nextobsf(&obss, &iobsu, 1)) <= 0) return -1;
		if (popt->intpref) {
			for (; (nr = nextobsf(&obss, &iobsr, 2)) > 0; iobsr += nr)
				if (timediff(obss.data[iobsr].time, obss.data[iobsu].time) > -DTTOL) break;
		}
		else {
			for (i = iobsr; (nr = nextobsf(&obss, &i, 2)) > 0; iobsr = i, i += nr)
				if (timediff(obss.data[i].time, obss.data[iobsu].time) > DTTOL) break;
		}
		nr = nextobsf(&obss, &iobsr, 2);
		if (nr <= 0) {
			nr = nextobsf(&obss, &iobsr, 2);
		}
		for (i = 0; i < nu && n < MAXOBS * 2; i++) obs[n++] = obss.data[iobsu + i]; /* 先存放rcv=1（移动站）的观测数据 */
		for (i = 0; i < nr && n < MAXOBS * 2; i++) obs[n++] = obss.data[iobsr + i]; /* 再存放rcv=2（基准站）的观测数据 */
		iobsu += nu;

		/* update sbas corrections */
		while (isbs < sbss.n) {
			time = gpst2time(sbss.msgs[isbs].week, sbss.msgs[isbs].tow);

			if (getbitu(sbss.msgs[isbs].msg, 8, 6) != 9) { /* except for geo nav */
				sbsupdatecorr(sbss.msgs + isbs, &navs);
			}
			if (timediff(time, obs[0].time) > -1.0 - DTTOL) break;
			isbs++;
		}
		/* update lex corrections */
		while (ilex < lexs.n) {
			if (lexupdatecorr(lexs.msgs + ilex, &navs, &time)) {
				if (timediff(time, obs[0].time) > -1.0 - DTTOL) break;
			}
			ilex++;
		}
		/* update rtcm ssr corrections */
		if (*rtcm_file) {
			update_rtcm_ssr(obs[0].time);
		}
	}
	else { /* input backward data */
		if ((nu = nextobsb(&obss, &iobsu, 1)) <= 0) return -1;
		if (popt->intpref) {
			for (; (nr = nextobsb(&obss, &iobsr, 2)) > 0; iobsr -= nr)
				if (timediff(obss.data[iobsr].time, obss.data[iobsu].time) < DTTOL) break;
		}
		else {
			for (i = iobsr; (nr = nextobsb(&obss, &i, 2)) > 0; iobsr = i, i -= nr)
				if (timediff(obss.data[i].time, obss.data[iobsu].time) < -DTTOL) break;
		}
		nr = nextobsb(&obss, &iobsr, 2);
		for (i = 0; i < nu && n < MAXOBS * 2; i++) obs[n++] = obss.data[iobsu - nu + 1 + i];
		for (i = 0; i < nr && n < MAXOBS * 2; i++) obs[n++] = obss.data[iobsr - nr + 1 + i];
		iobsu -= nu;

		/* update sbas corrections */
		while (isbs >= 0) {
			time = gpst2time(sbss.msgs[isbs].week, sbss.msgs[isbs].tow);

			if (getbitu(sbss.msgs[isbs].msg, 8, 6) != 9) { /* except for geo nav */
				sbsupdatecorr(sbss.msgs + isbs, &navs);
			}
			if (timediff(time, obs[0].time) < 1.0 + DTTOL) break;
			isbs--;
		}
		/* update lex corrections */
		while (ilex >= 0) {
			if (lexupdatecorr(lexs.msgs + ilex, &navs, &time)) {
				if (timediff(time, obs[0].time) < 1.0 + DTTOL) break;
			}
			ilex--;
		}
	}
	return n;
}
/* output to file message of invalid time mark -------------------------------*/
static void outinvalidtm(FILE* fptm, const solopt_t* opt, const gtime_t tm)
{
	gtime_t time = tm;
	double gpst;
	int week, timeu;
	char s[100];

	timeu = opt->timeu < 0 ? 0 : (opt->timeu > 20 ? 20 : opt->timeu);

	if (opt->times >= TIMES_UTC) time = gpst2utc(time);
	if (opt->times == TIMES_JST) time = timeadd(time, 9 * 3600.0);

	if (opt->timef) time2str(time, s, timeu);
	else {
		gpst = time2gpst(time, &week);
		if (86400 * 7 - gpst < 0.5 / pow(10.0, timeu)) {
			week++;
			gpst = 0.0;
		}
		sprintf(s, "%4d   %*.*f", week, 6 + (timeu <= 0 ? 0 : timeu + 1), timeu, gpst);
	}
	strcat(s, "   Q=0, Time mark is not valid\n");

	fwrite(s, strlen(s), 1, fptm);
}
/* fill structure sol_t for time mark ----------------------------------------*/
static sol_t fillsoltm(const sol_t solold, const sol_t solnew, const gtime_t tm)
{
	gtime_t t1 = { 0 }, t2 = { 0 };
	sol_t sol = solold;
	int i = 0;

	if (solold.stat == 0 || solnew.stat == 0) {
		sol.stat = 0;
	}
	else {
		sol.stat = (solold.stat > solnew.stat) ? solold.stat : solnew.stat;
	}
	sol.ns = (solold.ns < solnew.ns) ? solold.ns : solnew.ns;
	sol.ratio = (solold.ratio < solnew.ratio) ? solold.ratio : solnew.ratio;

	/* interpolation position and speed of time mark */
	t1 = solold.time;
	t2 = solnew.time;
	sol.time = tm;

	for (i = 0; i < 6; i++)
	{
		sol.rr[i] = solold.rr[i] + timediff(tm, t1) / timediff(t2, t1) * (solnew.rr[i] - solold.rr[i]);
	}

	return sol;
}
/* carrier-phase bias correction by fcb --------------------------------------*/
static void corr_phase_bias_fcb(obsd_t* obs, int n, const nav_t* nav)
{
	int i, j, k;

	for (i = 0; i < nav->nf; i++) {
		if (timediff(nav->fcb[i].te, obs[0].time) < -1E-3) continue;
		if (timediff(nav->fcb[i].ts, obs[0].time) > 1E-3) break;
		for (j = 0; j < n; j++) {
			for (k = 0; k < NFREQ; k++) {
				if (obs[j].L[k] == 0.0) continue;
				obs[j].L[k] -= nav->fcb[i].bias[obs[j].sat - 1][k];
			}
		}
		return;
	}
}
/* carrier-phase bias correction by ssr --------------------------------------*/
static void corr_phase_bias_ssr(obsd_t* obs, int n, const nav_t* nav)
{
	double lam;
	int i, j, code;

	for (i = 0; i < n; i++) for (j = 0; j < NFREQ; j++) {
		if (!(code = obs[i].code[j])) continue;
		if ((lam = nav->lam[obs[i].sat - 1][j]) == 0.0) continue;

		/* correct phase bias (cyc) */
		obs[i].L[j] -= nav->ssr[obs[i].sat - 1].pbias[code - 1] / lam;
	}
}
/* process positioning -------------------------------------------------------*/
static void procpos(FILE* fp, FILE* fptm, const prcopt_t* popt, const solopt_t* sopt, rtk_t* rtk, int mode)
{
	gtime_t time = { 0 };
	sol_t sol = { {0} }, oldsol = { {0} }, newsol = { {0} };
	obsd_t obs[MAXOBS * 2]; /* for rover and base */
	double rb[3] = { 0 };
	int i, nobs, n, solstatic, num = 0, pri[] = { 6,1,2,3,4,5,1,6 };

	trace(3, "procpos : mode=%d\n", mode);

	solstatic = sopt->solstatic &&
		(popt->mode == PMODE_STATIC || popt->mode == PMODE_STATIC_START || popt->mode == PMODE_PPP_STATIC);

	/* initialize unless running backwards on a combined run with continuous AR in which case keep the current states */
	if (mode == 0 || !revs || popt->modear == ARMODE_FIXHOLD)
		rtkinit(rtk, popt);

	rtcm_path[0] = '\0';

	/* 每次处理一个历元的数据，先输入观测数据，再计算站位置，最后输出结果 */
	while ((nobs = inputobs(obs, rtk->sol.stat, popt)) >= 0) {
		/* exclude satellites */
		for (i = n = 0; i < nobs; i++) {
			/* 判定当前观测值的卫星系统与后处理预设系统一致，且预设中的exsat选项不为1时，保留该观测值，否则剔除 */
			if ((satsys(obs[i].sat, NULL) & popt->navsys) && popt->exsats[obs[i].sat - 1] != 1) {
				obs[n++] = obs[i];
				//obs[i] = (obsd_t){ 0 }; // QHY添加
			}
		}
		if (n <= 0) continue;

		/* carrier-phase bias correction */
		if (navs.nf > 0) {
			corr_phase_bias_fcb(obs, n, &navs);
		}
		else if (!strstr(popt->pppopt, "-DIS_FCB")) {
			corr_phase_bias_ssr(obs, n, &navs);
		}
		/* disable L2 */
#if 0
		if (popt->freqopt == 1) {
			for (i = 0; i < n; i++) obs[i].L[1] = obs[i].P[1] = 0.0;
		}
#endif
		if (!rtkpos(rtk, obs, n, &navs)) {
			if (rtk->sol.eventime.time != 0) {
				if (mode == 0) {
					outinvalidtm(fptm, sopt, rtk->sol.eventime);
				}
				else if (!revs) {
					invalidtm[nitm++] = rtk->sol.eventime;
				}
			}
			continue;
		}

		if (mode == 0) { /* forward/backward */
			if (!solstatic) {
				outsol(fp, &rtk->sol, rtk->rb, sopt);
			}
			else if (time.time == 0 || pri[rtk->sol.stat] <= pri[sol.stat]) {
				sol = rtk->sol;
				for (i = 0; i < 3; i++) rb[i] = rtk->rb[i];
				if (time.time == 0 || timediff(rtk->sol.time, time) < 0.0) {
					time = rtk->sol.time;
				}
			}
			/* check time mark */
			if (rtk->sol.eventime.time != 0)
			{
				newsol = fillsoltm(oldsol, rtk->sol, rtk->sol.eventime);
				num++;
				if (!solstatic && mode == 0) {
					outsol(fptm, &newsol, rb, sopt);
				}
			}
			oldsol = rtk->sol;
		}
		else if (!revs) { /* combined-forward */
			if (isolf >= nepoch) return;
			solf[isolf] = rtk->sol;
			for (i = 0; i < 3; i++) rbf[i + isolf * 3] = rtk->rb[i];
			isolf++;
		}
		else { /* combined-backward */
			if (isolb >= nepoch) return;
			solb[isolb] = rtk->sol;
			for (i = 0; i < 3; i++) rbb[i + isolb * 3] = rtk->rb[i];
			isolb++;
		}
	}
	if (mode == 0 && solstatic && time.time != 0.0) {
		sol.time = time;
		outsol(fp, &sol, rb, sopt);
	}
}
/* validation of combined solutions ------------------------------------------*/
static int valcomb(const sol_t* solf, const sol_t* solb)
{
	double dr[3], var[3];
	int i;
	char tstr[32];

	trace(3, "valcomb :\n");

	/* compare forward and backward solution */
	for (i = 0; i < 3; i++) {
		dr[i] = solf->rr[i] - solb->rr[i];
		var[i] = solf->qr[i] + solb->qr[i];
	}
	for (i = 0; i < 3; i++) {
		if (dr[i] * dr[i] <= 16.0 * var[i]) continue; /* ok if in 4-sigma */

		time2str(solf->time, tstr, 2);
		trace(2, "degrade fix to float: %s dr=%.3f %.3f %.3f std=%.3f %.3f %.3f\n",
			tstr + 11, dr[0], dr[1], dr[2], SQRT(var[0]), SQRT(var[1]), SQRT(var[2]));
		return 0;
	}
	return 1;
}
/* combine forward/backward solutions and output results ---------------------*/
static void combres(FILE* fp, FILE* fptm, const prcopt_t* popt, const solopt_t* sopt)
{
	gtime_t time = { 0 };
	sol_t sols = { {0} }, sol = { {0} }, oldsol = { {0} }, newsol = { {0} };
	double tt, Qf[9], Qb[9], Qs[9], rbs[3] = { 0 }, rb[3] = { 0 }, rr_f[3], rr_b[3], rr_s[3];
	int i, j, k, solstatic, num = 0, pri[] = { 0,1,2,3,4,5,1,6 };

	trace(3, "combres : isolf=%d isolb=%d\n", isolf, isolb);

	solstatic = sopt->solstatic &&
		(popt->mode == PMODE_STATIC || popt->mode == PMODE_STATIC_START || popt->mode == PMODE_PPP_STATIC);

	for (i = 0, j = isolb - 1; i < isolf && j >= 0; i++, j--) {
		if ((tt = timediff(solf[i].time, solb[j].time)) < -DTTOL) {
			sols = solf[i];
			for (k = 0; k < 3; k++) rbs[k] = rbf[k + i * 3];
			j++;
		}
		else if (tt > DTTOL) {
			sols = solb[j];
			for (k = 0; k < 3; k++) rbs[k] = rbb[k + j * 3];
			i--;
		}
		else if (solf[i].stat < solb[j].stat) {
			sols = solf[i];
			for (k = 0; k < 3; k++) rbs[k] = rbf[k + i * 3];
		}
		else if (solf[i].stat > solb[j].stat) {
			sols = solb[j];
			for (k = 0; k < 3; k++) rbs[k] = rbb[k + j * 3];
		}
		else {
			sols = solf[i];
			sols.time = timeadd(sols.time, -tt / 2.0);

			if ((popt->mode == PMODE_KINEMA || popt->mode == PMODE_MOVEB) &&
				sols.stat == SOLQ_FIX) {
				/* degrade fix to float if validation failed */
				if (!valcomb(solf + i, solb + j)) sols.stat = SOLQ_FLOAT;
			}
			for (k = 0; k < 3; k++) {
				Qf[k + k * 3] = solf[i].qr[k];
				Qb[k + k * 3] = solb[j].qr[k];
			}
			Qf[1] = Qf[3] = solf[i].qr[3];
			Qf[5] = Qf[7] = solf[i].qr[4];
			Qf[2] = Qf[6] = solf[i].qr[5];
			Qb[1] = Qb[3] = solb[j].qr[3];
			Qb[5] = Qb[7] = solb[j].qr[4];
			Qb[2] = Qb[6] = solb[j].qr[5];

			if (popt->mode == PMODE_MOVEB) {
				for (k = 0; k < 3; k++) rr_f[k] = solf[i].rr[k] - rbf[k + i * 3];
				for (k = 0; k < 3; k++) rr_b[k] = solb[j].rr[k] - rbb[k + j * 3];
				if (smoother(rr_f, Qf, rr_b, Qb, 3, rr_s, Qs)) continue;
				for (k = 0; k < 3; k++) sols.rr[k] = rbs[k] + rr_s[k];
			}
			else {
				if (smoother(solf[i].rr, Qf, solb[j].rr, Qb, 3, sols.rr, Qs)) continue;
			}
			sols.qr[0] = (float)Qs[0];
			sols.qr[1] = (float)Qs[4];
			sols.qr[2] = (float)Qs[8];
			sols.qr[3] = (float)Qs[1];
			sols.qr[4] = (float)Qs[5];
			sols.qr[5] = (float)Qs[2];

			/* smoother for velocity solution */
			if (popt->dynamics) {
				for (k = 0; k < 3; k++) {
					Qf[k + k * 3] = solf[i].qv[k];
					Qb[k + k * 3] = solb[j].qv[k];
				}
				Qf[1] = Qf[3] = solf[i].qv[3];
				Qf[5] = Qf[7] = solf[i].qv[4];
				Qf[2] = Qf[6] = solf[i].qv[5];
				Qb[1] = Qb[3] = solb[j].qv[3];
				Qb[5] = Qb[7] = solb[j].qv[4];
				Qb[2] = Qb[6] = solb[j].qv[5];
				if (smoother(solf[i].rr + 3, Qf, solb[j].rr + 3, Qb, 3, sols.rr + 3, Qs)) continue;
				sols.qv[0] = (float)Qs[0];
				sols.qv[1] = (float)Qs[4];
				sols.qv[2] = (float)Qs[8];
				sols.qv[3] = (float)Qs[1];
				sols.qv[4] = (float)Qs[5];
				sols.qv[5] = (float)Qs[2];
			}
		}
		if (!solstatic) {
			outsol(fp, &sols, rbs, sopt);
		}
		else if (time.time == 0 || pri[sols.stat] <= pri[sol.stat]) {
			sol = sols;
			for (k = 0; k < 3; k++) rb[k] = rbs[k];
			if (time.time == 0 || timediff(sols.time, time) < 0.0) {
				time = sols.time;
			}
		}
		if (iitm < nitm && timediff(invalidtm[iitm], sols.time) < 0.0)
		{
			outinvalidtm(fptm, sopt, invalidtm[iitm]);
			iitm++;
		}
		if (sols.eventime.time != 0)
		{
			newsol = fillsoltm(oldsol, sols, sols.eventime);
			num++;
			if (!solstatic) {
				outsol(fptm, &newsol, rb, sopt);
			}
		}
		oldsol = sols;
	}
	if (solstatic && time.time != 0.0) {
		sol.time = time;
		outsol(fp, &sol, rb, sopt);
	}
}
/* read prec ephemeris, sbas data, lex data, tec grid and open rtcm ----------*/
static void readpreceph(char** infile, int n, const prcopt_t* prcopt,
	nav_t* nav, sbs_t* sbs, lex_t* lex)
{
	seph_t seph0 = { 0 };
	int i;
	char* ext;

	trace(2, "readpreceph: n=%d\n", n);

	nav->ne = nav->nemax = 0;
	nav->nc = nav->ncmax = 0;
	nav->nf = nav->nfmax = 0;
	sbs->n = sbs->nmax = 0;
	lex->n = lex->nmax = 0;

	/* read precise ephemeris files */
	for (i = 0; i < n; i++) {
		if (strstr(infile[i], "%r") || strstr(infile[i], "%b")) continue;
		readsp3(infile[i], nav, 0);
	}
	/* read precise clock files */
	for (i = 0; i < n; i++) {
		if (strstr(infile[i], "%r") || strstr(infile[i], "%b")) continue;
		readrnxc(infile[i], nav);
	}
	/* read satellite fcb files */
	for (i = 0; i < n; i++) {
		if (strstr(infile[i], "%r") || strstr(infile[i], "%b")) continue;
		if ((ext = strrchr(infile[i], '.')) &&
			(!strcmp(ext, ".fcb") || !strcmp(ext, ".FCB"))) {
			readfcb(infile[i], nav);
		}
	}
	/* read solution status files for ppp correction */
	for (i = 0; i < n; i++) {
		if (strstr(infile[i], "%r") || strstr(infile[i], "%b")) continue;
		if ((ext = strrchr(infile[i], '.')) &&
			(!strcmp(ext, ".stat") || !strcmp(ext, ".STAT") ||
				!strcmp(ext, ".stec") || !strcmp(ext, ".STEC") ||
				!strcmp(ext, ".trp") || !strcmp(ext, ".TRP"))) {
			pppcorr_read(&nav->pppcorr, infile[i]);
		}
	}
	/* read sbas message files */
	for (i = 0; i < n; i++) {
		if (strstr(infile[i], "%r") || strstr(infile[i], "%b")) continue;
		sbsreadmsg(infile[i], prcopt->sbassatsel, sbs);
	}
	/* read lex message files */
	for (i = 0; i < n; i++) {
		if (strstr(infile[i], "%r") || strstr(infile[i], "%b")) continue;
		lexreadmsg(infile[i], 0, lex);
	}
	/* allocate sbas ephemeris */
	nav->ns = nav->nsmax = NSATSBS * 2;
	if (!(nav->seph = (seph_t*)malloc(sizeof(seph_t) * nav->ns))) {
		showmsg("error : sbas ephem memory allocation");
		trace(1, "error : sbas ephem memory allocation");
		return;
	}
	for (i = 0; i < nav->ns; i++) nav->seph[i] = seph0;

	/* set rtcm file and initialize rtcm struct */
	rtcm_file[0] = rtcm_path[0] = '\0'; fp_rtcm = NULL;

	for (i = 0; i < n; i++) {
		if ((ext = strrchr(infile[i], '.')) &&
			(!strcmp(ext, ".rtcm3") || !strcmp(ext, ".RTCM3"))) {
			strcpy(rtcm_file, infile[i]);
			init_rtcm(&rtcm);
			break;
		}
	}
}
/* free prec ephemeris and sbas data -----------------------------------------*/
static void freepreceph(nav_t* nav, sbs_t* sbs, lex_t* lex)
{
	int i;

	trace(3, "freepreceph:\n");

	free(nav->peph); nav->peph = NULL; nav->ne = nav->nemax = 0;
	free(nav->pclk); nav->pclk = NULL; nav->nc = nav->ncmax = 0;
	free(nav->fcb); nav->fcb = NULL; nav->nf = nav->nfmax = 0;
	free(nav->seph); nav->seph = NULL; nav->ns = nav->nsmax = 0;
	free(sbs->msgs); sbs->msgs = NULL; sbs->n = sbs->nmax = 0;
	free(lex->msgs); lex->msgs = NULL; lex->n = lex->nmax = 0;
	for (i = 0; i < nav->nt; i++) {
		free(nav->tec[i].data);
		free(nav->tec[i].rms);
	}
	free(nav->tec); nav->tec = NULL; nav->nt = nav->ntmax = 0;

	if (fp_rtcm) fclose(fp_rtcm);
	free_rtcm(&rtcm);
}
/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, char** infile,
	const int* index, int n, const prcopt_t* prcopt,
	obs_t* obs, nav_t* nav, sta_t* sta)
{
	int i, j, ind = 0, nobs = 0, rcv = 1;

	trace(3, "readobsnav: ts=%s n=%d\n", time_str(ts, 0), n);

	/* 初始化obs,nav结构体中的参数 */
	obs->data = NULL; obs->n = obs->nmax = 0;
	nav->eph = NULL; nav->n = nav->nmax = 0;
	nav->geph = NULL; nav->ng = nav->ngmax = 0;
	/* free(nav->seph); */ /* is this needed to avoid memory leak??? */
	nav->seph = NULL; nav->ns = nav->nsmax = 0;
	nepoch = 0; /* number of observation epochs */

	/* 通过for循环遍历次数n，读取rnx文件数据到nav，obs中 */
	for (i = 0; i < n; i++) {
		if (checkbrk("")) return 0;

		if (index[i] != ind) {
			if (obs->n > nobs) rcv++;
			ind = index[i]; nobs = obs->n; /* 通过index,obs给当前ind,nobs赋值*/
		}
		/* read rinex obs and nav file */
		if (readrnxt(infile[i], rcv, ts, te, ti, prcopt->rnxopt[rcv <= 1 ? 0 : 1], obs, nav,
			rcv <= 2 ? sta + rcv - 1 : NULL) < 0) {
			checkbrk("error : insufficient memory");
			trace(1, "insufficient memory\n");
			return 0;
		}
	}
	if (obs->n <= 0) {
		checkbrk("error : no obs data");
		trace(1, "\n");
		return 0;
	}
	if (nav->n <= 0 && nav->ng <= 0 && nav->ns <= 0) {
		checkbrk("error : no nav data");
		trace(1, "\n");
		return 0;
	}
	/* sort observation data 给obs排序 */
	nepoch = sortobs(obs);

	/* delete duplicated ephemeris */
	uniqnav(nav);

	/* Set lam[1] as flag for satellite PCV antenna offset calc to indicate
	   using Galileo E5a or E5b  */
	if (prcopt->nf > 2) {       /* if L1+L2+L5 solution */
		for (i = 0; i < MAXSAT; i++) {
			/* set Galileo lam[1]=0 which will force offset calc to use lam[2] (E5a) */
			if (satsys(i + 1, NULL) == SYS_GAL) nav->lam[i][1] = 0;
		}
	}

	/* set time span for progress display */
	if (ts.time == 0 || te.time == 0) {
		for (i = 0; i < obs->n; i++) if (obs->data[i].rcv == 1) break; /* 搜寻接收机号rcv为1的第一个obs位置索引i */
		for (j = obs->n - 1; j >= 0; j--) if (obs->data[j].rcv == 1) break; /* 搜寻接收机号rcv为1的最后一个obs位置索引j */
		if (i < j) {
			if (ts.time == 0) ts = obs->data[i].time;
			if (te.time == 0) te = obs->data[j].time;
			settspan(ts, te);
		}
	}
	return 1;
}
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t* obs, nav_t* nav)
{
	trace(3, "freeobsnav:\n");

	free(obs->data); obs->data = NULL; obs->n = obs->nmax = 0;
	free(nav->eph); nav->eph = NULL; nav->n = nav->nmax = 0;
	free(nav->geph); nav->geph = NULL; nav->ng = nav->ngmax = 0;
	free(nav->seph); nav->seph = NULL; nav->ns = nav->nsmax = 0;
}
/* average of single position ------------------------------------------------*/
static int avepos(double* ra, int rcv, const obs_t* obs, const nav_t* nav,
	const prcopt_t* opt)
{
	obsd_t data[MAXOBS];
	gtime_t ts = { 0 };
	sol_t sol = { {0} };
	int i, j, n = 0, m, iobs;
	char msg[128];

	trace(3, "avepos: rcv=%d obs.n=%d\n", rcv, obs->n);

	for (i = 0; i < 3; i++) ra[i] = 0.0; /* 先给ra赋初值0 */

	for (iobs = 0; (m = nextobsf(obs, &iobs, rcv)) > 0; iobs += m) {
		for (i = j = 0; i < m && i < MAXOBS; i++) {
			data[j] = obs->data[iobs + i];
			if ((satsys(data[j].sat, NULL) & opt->navsys) && opt->exsats[data[j].sat - 1] != 1) j++;
		}
		if (j <= 0 || !screent(data[0].time, ts, ts, 1.0)) continue; /* only 1 hz */

		if (!pntpos(data, j, nav, opt, &sol, NULL, NULL, msg)) continue;

		for (i = 0; i < 3; i++) ra[i] += sol.rr[i];
		n++;
	}
	if (n <= 0) {
		trace(1, "no average of base station position\n");
		return 0;
	}
	for (i = 0; i < 3; i++) ra[i] /= n;
	return 1;
}
/* station position from file ------------------------------------------------*/
static int getstapos(const char* file, char* name, double* r)
{
	FILE* fp;
	char buff[256], sname[256], * p, * q;
	double pos[3];

	trace(3, "getstapos: file=%s name=%s\n", file, name);

	if (!(fp = fopen(file, "r"))) {
		trace(1, "station position file open error: %s\n", file);
		return 0;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		if ((p = strchr(buff, '%'))) *p = '\0';

		if (sscanf(buff, "%lf %lf %lf %s", pos, pos + 1, pos + 2, sname) < 4) continue;

		for (p = sname, q = name; *p && *q; p++, q++) {
			if (toupper((int)*p) != toupper((int)*q)) break;
		}
		if (!*p) {
			pos[0] *= D2R;
			pos[1] *= D2R;
			pos2ecef(pos, r);
			fclose(fp);
			return 1;
		}
	}
	fclose(fp);
	trace(1, "no station position: %s %s\n", name, file);
	return 0;
}
/* antenna phase center position ---------------------------------------------*/
static int antpos(prcopt_t* opt, int rcvno, const obs_t* obs, const nav_t* nav,
	const sta_t* sta, const char* posfile)
{
	double* rr = rcvno == 1 ? opt->ru : opt->rb, del[3], pos[3], dr[3] = { 0 };
	int i, postype = rcvno == 1 ? opt->rovpos : opt->refpos;
	char* name;

	trace(3, "antpos  : rcvno=%d\n", rcvno);

	if (postype == POSOPT_SINGLE) { /* average of single position */
		if (!avepos(rr, rcvno, obs, nav, opt)) {
			showmsg("error : station pos computation");
			return 0;
		}
	}
	else if (postype == POSOPT_FILE) { /* read from position file */
		name = stas[rcvno == 1 ? 0 : 1].name;
		if (!getstapos(posfile, name, rr)) {
			showmsg("error : no position of %s in %s", name, posfile);
			return 0;
		}
	}
	else if (postype == POSOPT_RINEX) { /* get from rinex header */
		if (norm(stas[rcvno == 1 ? 0 : 1].pos, 3) <= 0.0) {
			showmsg("error : no position in rinex header");
			trace(1, "no position in rinex header\n");
			return 0;
		}
		/* add antenna delta unless already done in antpcv() */
		if (!strcmp(opt->anttype[rcvno], "*")) {
			if (stas[rcvno == 1 ? 0 : 1].deltype == 0) { /* enu */
				for (i = 0; i < 3; i++) del[i] = stas[rcvno == 1 ? 0 : 1].del[i];
				del[2] += stas[rcvno == 1 ? 0 : 1].hgt;
				ecef2pos(stas[rcvno == 1 ? 0 : 1].pos, pos);
				enu2ecef(pos, del, dr);
			}
			else { /* xyz */
				for (i = 0; i < 3; i++) dr[i] = stas[rcvno == 1 ? 0 : 1].del[i];
			}
		}
		for (i = 0; i < 3; i++) rr[i] = stas[rcvno == 1 ? 0 : 1].pos[i] + dr[i];
	}
	return 1;
}
/* open procssing session ----------------------------------------------------*/
static int openses(const prcopt_t* popt, const solopt_t* sopt,
	const filopt_t* fopt, nav_t* nav, pcvs_t* pcvs, pcvs_t* pcvr)
{
	int i;

	trace(3, "openses :\n");

	/* read satellite antenna parameters */
	if (*fopt->satantp && !(readpcv(fopt->satantp, pcvs))) {
		showmsg("error : no sat ant pcv in %s", fopt->satantp);
		trace(1, "sat antenna pcv read error: %s\n", fopt->satantp);
		return 0;
	}
	/* read receiver antenna parameters */
	if (*fopt->rcvantp && !(readpcv(fopt->rcvantp, pcvr))) {
		showmsg("error : no rec ant pcv in %s", fopt->rcvantp);
		trace(1, "rec antenna pcv read error: %s\n", fopt->rcvantp);
		return 0;
	}
	/* open geoid data */
	if (sopt->geoid > 0 && *fopt->geoid) {
		if (!opengeoid(sopt->geoid, fopt->geoid)) {
			showmsg("error : no geoid data %s", fopt->geoid);
			trace(2, "no geoid data %s\n", fopt->geoid);
		}
	}
	/* use satellite L2 offset if L5 offset does not exists */
	for (i = 0; i < pcvs->n; i++) {
		if (norm(pcvs->pcv[i].off[2], 3) > 0.0) continue;
		matcpy(pcvs->pcv[i].off[2], pcvs->pcv[i].off[1], 3, 1);
		matcpy(pcvs->pcv[i].var[2], pcvs->pcv[i].var[1], 19, 1);
	}
	for (i = 0; i < pcvr->n; i++) {
		if (norm(pcvr->pcv[i].off[2], 3) > 0.0) continue;
		matcpy(pcvr->pcv[i].off[2], pcvr->pcv[i].off[1], 3, 1);
		matcpy(pcvr->pcv[i].var[2], pcvr->pcv[i].var[1], 19, 1);
	}
	return 1;
}
/* close procssing session ---------------------------------------------------*/
static void closeses(nav_t* nav, pcvs_t* pcvs, pcvs_t* pcvr)
{
	trace(3, "closeses:\n");

	/* free antenna parameters */
	free(pcvs->pcv); pcvs->pcv = NULL; pcvs->n = pcvs->nmax = 0;
	free(pcvr->pcv); pcvr->pcv = NULL; pcvr->n = pcvr->nmax = 0;

	/* close geoid data */
	closegeoid();

	/* free erp data */
	free(nav->erp.data); nav->erp.data = NULL; nav->erp.n = nav->erp.nmax = 0;

	/* close solution statistics and debug trace */
	rtkclosestat();
	traceclose();
}
/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t* popt, nav_t* nav, const pcvs_t* pcvs,
	const pcvs_t* pcvr, const sta_t* sta)
{
	pcv_t* pcv;
	double pos[3], del[3];
	int i, j, mode = PMODE_DGPS <= popt->mode && popt->mode <= PMODE_FIXED;
	char id[64];

	/* set satellite antenna parameters */
	for (i = 0; i < MAXSAT; i++) {
		if (!(satsys(i + 1, NULL) & popt->navsys)) continue;
		if (!(pcv = searchpcv(i + 1, "", time, pcvs))) {
			satno2id(i + 1, id);
			trace(3, "no satellite antenna pcv: %s\n", id);
			continue;
		}
		nav->pcvs[i] = *pcv;
	}
	for (i = 0; i < (mode ? 2 : 1); i++) {
		if (!strcmp(popt->anttype[i], "*")) { /* set by station parameters */
			strcpy(popt->anttype[i], sta[i].antdes);
			if (sta[i].deltype == 1) { /* xyz */
				if (norm(sta[i].pos, 3) > 0.0) {
					ecef2pos(sta[i].pos, pos);
					ecef2enu(pos, sta[i].del, del);
					for (j = 0; j < 3; j++) popt->antdel[i][j] = del[j];
				}
			}
			else { /* enu */
				for (j = 0; j < 3; j++) popt->antdel[i][j] = stas[i].del[j];
			}
		}
		if (!(pcv = searchpcv(0, popt->anttype[i], time, pcvr))) {
			trace(2, "no receiver antenna pcv: %s\n", popt->anttype[i]);
			*popt->anttype[i] = '\0';
			continue;
		}
		strcpy(popt->anttype[i], pcv->type);
		popt->pcvr[i] = *pcv;
	}
}
/* read ocean tide loading parameters ----------------------------------------*/
static void readotl(prcopt_t* popt, const char* file, const sta_t* sta)
{
	int i, mode = PMODE_DGPS <= popt->mode && popt->mode <= PMODE_FIXED;

	for (i = 0; i < (mode ? 2 : 1); i++) {
		readblq(file, sta[i].name, popt->odisp[i]);
	}
}
/* write header to output file -----------------------------------------------*/
static int outhead(const char* outfile, char** infile, int n,
	const prcopt_t* popt, const solopt_t* sopt)
{
	FILE* fp = stdout;

	trace(3, "outhead: outfile=%s n=%d\n", outfile, n);

	if (*outfile) {
		createdir(outfile);

		if (!(fp = fopen(outfile, "w"))) {
			showmsg("error : open output file %s", outfile);
			return 0;
		}
	}
	/* output header */
	outheader(fp, infile, n, popt, sopt);

	if (*outfile) fclose(fp);

	return 1;
}
/* open output file for append -----------------------------------------------*/
static FILE* openfile(const char* outfile)
{
	trace(3, "openfile: outfile=%s\n", outfile);

	return !*outfile ? stdout : fopen(outfile, "a");
}
/* Name time marks file ------------------------------------------------------*/
static void namefiletm(char* outfiletm, const char* outfile)
{
	int i;

	for (i = strlen(outfile); i > 0; i--) {
		if (outfile[i] == '.') {
			break;
		}
	}
	/* if no file extension, then name time marks file as name of outfile + _events.pos */
	if (i == 0) {
		i = strlen(outfile);
	}
	strncpy(outfiletm, outfile, i);
	strcat(outfiletm, "_events.pos");
}
/* execute processing session ------------------------------------------------*/
static int execses(gtime_t ts, gtime_t te, double ti, const prcopt_t* popt,
	const solopt_t* sopt, const filopt_t* fopt, int flag,
	char** infile, const int* index, int n, char* outfile)
{
	FILE* fp, * fptm;
	rtk_t rtk;
	prcopt_t popt_ = *popt;
	solopt_t tmsopt = *sopt;
	char tracefile[1024], statfile[1024], path[1024], * ext, outfiletm[1024] = { 0 };
	int i, j, k;

	trace(3, "execses : n=%d outfile=%s\n", n, outfile);

	/* open debug trace */
	if (flag && sopt->trace > 0) {
		if (*outfile) {
			strcpy(tracefile, outfile);
			strcat(tracefile, ".trace");
		}
		else {
			strcpy(tracefile, fopt->trace);
		}
		traceclose();
		traceopen(tracefile);
		tracelevel(sopt->trace);
	}
	/* read ionosphere data file */
	if (*fopt->iono && (ext = strrchr(fopt->iono, '.'))) {
		if (strlen(ext) == 4 && (ext[3] == 'i' || ext[3] == 'I')) {
			reppath(fopt->iono, path, ts, "", "");
			readtec(path, &navs, 1);
		}
	}
	/* read erp data */
	if (*fopt->eop) {
		free(navs.erp.data); navs.erp.data = NULL; navs.erp.n = navs.erp.nmax = 0;
		reppath(fopt->eop, path, ts, "", "");
		if (!readerp(path, &navs.erp)) {
			showmsg("error : no erp data %s", path);
			trace(2, "no erp data %s\n", path);
		}
	}
	/* read obs and nav data */
	if (!readobsnav(ts, te, ti, infile, index, n, &popt_, &obss, &navs, stas)) {
		/* free obs and nav data */
		freeobsnav(&obss, &navs);
		return 0;
	}
	/* QHY Functions ------------------------------------------ */
	eph_num = caleph(&obss);
	/*sele_satsys(&obss, &popt_, &stec_list);
	detslp(&stec_list);
	calSTEC(&stec_list);*/
	/* -------------------------------------------------------- */
	/* read dcb parameters */
	if (*fopt->dcb) {
		reppath(fopt->dcb, path, ts, "", "");
		readdcb(path, &navs, stas);
	}
	else {
		for (i = 0; i < 3; i++) {
			for (j = 0; j < MAXSAT; j++) navs.cbias[j][i] = 0;
			for (j = 0; j < MAXRCV; j++) for (k = 0; k < 2; k++) navs.rbias[j][k][i] = 0;
		}
	}
	/* set antenna parameters */
	if (popt_.mode != PMODE_SINGLE) {
		setpcv(obss.n > 0 ? obss.data[0].time : timeget(), &popt_, &navs, &pcvss, &pcvsr,
			stas);
	}
	/* read ocean tide loading parameters */
	if (popt_.mode > PMODE_SINGLE && *fopt->blq) {
		readotl(&popt_, fopt->blq, stas);
	}
	/* rover/reference fixed position */
	if (popt_.mode == PMODE_FIXED) {
		if (!antpos(&popt_, 1, &obss, &navs, stas, fopt->stapos)) {
			freeobsnav(&obss, &navs);
			return 0;
		}
		if (!antpos(&popt_, 2, &obss, &navs, stas, fopt->stapos)) {
			freeobsnav(&obss, &navs);
			return 0;
		}
	}
	else if (PMODE_DGPS <= popt_.mode && popt_.mode <= PMODE_STATIC_START) { /* RTK动态/静态定位基站位置计算入口 */
		if (!antpos(&popt_, 2, &obss, &navs, stas, fopt->stapos)) {
			freeobsnav(&obss, &navs);
			return 0;
		}
	}
	/* open solution statistics */
	if (flag && sopt->sstat > 0) {
		strcpy(statfile, outfile);
		strcat(statfile, ".stat");
		rtkclosestat();
		rtkopenstat(statfile, sopt->sstat);
	}
	/* write header to output file */
	if (flag && !outhead(outfile, infile, n, &popt_, sopt)) {
		freeobsnav(&obss, &navs);
		return 0;
	}
	/* name time events file */
	namefiletm(outfiletm, outfile);
	/* write header to file with time marks */
	outhead(outfiletm, infile, n, &popt_, &tmsopt);

	iobsu = iobsr = isbs = ilex = revs = aborts = 0;

	if (popt_.mode == PMODE_SINGLE || popt_.soltype == 0) {
		if ((fp = openfile(outfile)) && (fptm = openfile(outfiletm))) {
			procpos(fp, fptm, &popt_, sopt, &rtk, 0); /* forward */
			fclose(fp);
			fclose(fptm);
		}
	}
	else if (popt_.soltype == 1) {
		if ((fp = openfile(outfile)) && (fptm = openfile(outfiletm))) {
			revs = 1; iobsu = iobsr = obss.n - 1; isbs = sbss.n - 1; ilex = lexs.n - 1;
			procpos(fp, fptm, &popt_, sopt, &rtk, 0); /* backward */
			fclose(fp);
			fclose(fptm);
		}
	}
	else { /* combined */
		solf = (sol_t*)malloc(sizeof(sol_t) * nepoch);
		solb = (sol_t*)malloc(sizeof(sol_t) * nepoch);
		rbf = (double*)malloc(sizeof(double) * nepoch * 3);
		rbb = (double*)malloc(sizeof(double) * nepoch * 3);

		if (solf && solb) {
			isolf = isolb = 0;
			procpos(NULL, NULL, &popt_, sopt, &rtk, 1); /* forward */
			revs = 1; iobsu = iobsr = obss.n - 1; isbs = sbss.n - 1; ilex = lexs.n - 1;
			procpos(NULL, NULL, &popt_, sopt, &rtk, 1); /* backward */

			/* combine forward/backward solutions */
			if (!aborts && (fp = openfile(outfile)) && (fptm = openfile(outfiletm))) {
				combres(fp, fptm, &popt_, sopt);
				fclose(fp);
				fclose(fptm);
			}
		}
		else showmsg("error : memory allocation");
		free(solf);
		free(solb);
		free(rbf);
		free(rbb);
	}
	/* free rtk, obs and nav data */
	rtkfree(&rtk);
	freeobsnav(&obss, &navs);

	return aborts ? 1 : 0;
}
/* execute processing session for each rover ---------------------------------*/
static int execses_r(gtime_t ts, gtime_t te, double ti, const prcopt_t* popt,
	const solopt_t* sopt, const filopt_t* fopt, int flag,
	char** infile, const int* index, int n, char* outfile,
	const char* rov)
{
	gtime_t t0 = { 0 };
	int i, stat = 0;
	char* ifile[MAXINFILE], ofile[1024], * rov_, * p, * q, s[64] = "";

	trace(3, "execses_r: n=%d outfile=%s\n", n, outfile);

	for (i = 0; i < n; i++) if (strstr(infile[i], "%r")) break;

	if (i < n) { /* include rover keywords */
		if (!(rov_ = (char*)malloc(strlen(rov) + 1))) return 0;
		strcpy(rov_, rov);

		for (i = 0; i < n; i++) {
			if (!(ifile[i] = (char*)malloc(1024))) {
				free(rov_); for (; i >= 0; i--) free(ifile[i]);
				return 0;
			}
		}
		for (p = rov_;; p = q + 1) { /* for each rover */
			if ((q = strchr(p, ' '))) *q = '\0';

			if (*p) {
				strcpy(proc_rov, p);
				if (ts.time) time2str(ts, s, 0); else *s = '\0';
				if (checkbrk("reading    : %s", s)) {
					stat = 1;
					break;
				}
				for (i = 0; i < n; i++) reppath(infile[i], ifile[i], t0, p, "");
				reppath(outfile, ofile, t0, p, "");

				/* execute processing session */
				stat = execses(ts, te, ti, popt, sopt, fopt, flag, ifile, index, n, ofile);
			}
			if (stat == 1 || !q) break;
		}
		free(rov_); for (i = 0; i < n; i++) free(ifile[i]);
	}
	else {
		/* execute processing session */
		stat = execses(ts, te, ti, popt, sopt, fopt, flag, infile, index, n, outfile);
	}
	return stat;
}
/* execute processing session for each base station --------------------------*/
static int execses_b(gtime_t ts, gtime_t te, double ti, const prcopt_t* popt,
	const solopt_t* sopt, const filopt_t* fopt, int flag,
	char** infile, const int* index, int n, char* outfile,
	const char* rov, const char* base)
{
	gtime_t t0 = { 0 };
	int i, stat = 0;
	char* ifile[MAXINFILE], ofile[1024], * base_, * p, * q, s[64];

	trace(3, "execses_b: n=%d outfile=%s\n", n, outfile);

	/* read prec ephemeris and sbas data */
	readpreceph(infile, n, popt, &navs, &sbss, &lexs);

	for (i = 0; i < n; i++) if (strstr(infile[i], "%b")) break;

	if (i < n) { /* include base station keywords */
		if (!(base_ = (char*)malloc(strlen(base) + 1))) { /* 给base_分配内存地址 */
			freepreceph(&navs, &sbss, &lexs);
			return 0;
		}
		strcpy(base_, base);

		for (i = 0; i < n; i++) {
			if (!(ifile[i] = (char*)malloc(1024))) { /* 给每个ifile分配内存地址 */
				free(base_); for (; i >= 0; i--) free(ifile[i]);
				freepreceph(&navs, &sbss, &lexs);
				return 0;
			}
		}
		for (p = base_;; p = q + 1) { /* for each base station */
			if ((q = strchr(p, ' '))) *q = '\0'; /* q指向base中的空格位置，并将其转变为结束符 */

			if (*p) {
				strcpy(proc_base, p);
				if (ts.time) time2str(ts, s, 0); else *s = '\0';
				if (checkbrk("reading    : %s", s)) {
					stat = 1;
					break;
				}
				for (i = 0; i < n; i++) reppath(infile[i], ifile[i], t0, "", p);
				reppath(outfile, ofile, t0, "", p);

				stat = execses_r(ts, te, ti, popt, sopt, fopt, flag, ifile, index, n, ofile, rov);
			}
			if (stat == 1 || !q) break;
		}
		free(base_); for (i = 0; i < n; i++) free(ifile[i]);
	}
	else {
		stat = execses_r(ts, te, ti, popt, sopt, fopt, flag, infile, index, n, outfile, rov);
	}
	/* free prec ephemeris and sbas data */
	freepreceph(&navs, &sbss, &lexs);

	return stat;
}
/* QHY Defined Function ------------------------------------------------------*/
static int caleph(obs_t* obs) { /* 计算obss中历元个数并返回 */
	int i, j;
	for (i = 0, j = 0; i < obs->n; i++) {
		if (obs->data[i].time.time != obs->data[i + 1].time.time)j++;
	}
	return j;
}

static double cal_MW(double L1, double L2, double P1, double P2) {
	if (L1 == 0 || L2 == 0 || P1 == 0 || P2 == 0)return 0.0;
	double L_wl = a_wl * lam_1 * L1 - b_wl * lam_2 * L2;
	double P_nl = a_nl * P1 + b_nl * P2;
	double N_mw = (L_wl - P_nl) / lam_wl;
	return N_mw;
}

static double cal_GF(double L1, double L2, double P1, double P2) {
	if (L1 == 0 || L2 == 0 || P1 == 0 || P2 == 0)return 0.0;
	double L_gf = lam_1 * L1 - lam_2 * L2;
	double N_gf = L_gf / lam_1;
	return N_gf;
}

static void sele_satsys(obs_t* obs, prcopt_t* popt, stec_q* stecl) {
	obsd_t obs_temp[NSATCMP];
	malloc(stecl->stec, sizeof(stecd_q*) * eph_num);

	//system("pause");

	for (int i = 0; i < eph_num; i++) {
		int nu = 0; int n = 0; int n2 = 0;
		if ((nu = nextobsf(obs, &obs_idx, 1)) <= 0) return;
		for (int j = 0; j < nu && n < NSATCMP; j++) obs_temp[n++] = obs->data[obs_idx + j];
		stecl->stec[i].time = obs->data[obs_idx].time;
		obs_idx += nu;
		for (int k = 0; k < n; k++) {
			if ((satsys(obs_temp[k].sat, NULL) & popt->navsys) && (int)popt->exsats[obs_temp[k].sat - 1] != 1) {
				obs_temp[n2++] = obs_temp[k];
				//obs_temp[k] = (obsd_t){ 0 };
			}
		}
		for (int p = 0; p < NSATCMP; p++) {
			int flag = 1; /* 卫星数据存在标记符 */
			int brk = 0; if (i == 282 && p == 5)brk++;
			for (int q = 0; q < n2; q++) {
				if ((obs_temp[q].sat - SATCMP) == p) { /* 如果当前历元存在该卫星数据，则传递参数 */
					stecl->stec[i].L1[p] = obs_temp[q].L[0];
					stecl->stec[i].L2[p] = obs_temp[q].L[1];
					stecl->stec[i].P1[p] = obs_temp[q].P[0];
					stecl->stec[i].P2[p] = obs_temp[q].P[1];
					stecl->stec[i].N_mw[p] = cal_MW(stecl->stec[i].L1[p], stecl->stec[i].L2[p], stecl->stec[i].P1[p], stecl->stec[i].P2[p]);
					stecl->stec[i].N_gf[p] = cal_GF(stecl->stec[i].L1[p], stecl->stec[i].L2[p], stecl->stec[i].P1[p], stecl->stec[i].P2[p]);
					flag = 0; break;
				}
			}
			if (flag) { /* 无该卫星参数需将数据置0 */
				stecl->stec[i].L1[p] = 0.0;
				stecl->stec[i].L2[p] = 0.0;
				stecl->stec[i].P1[p] = 0.0;
				stecl->stec[i].P2[p] = 0.0;
			}
		}
		//printf("L1 = %f\n", (stec_tmp + i)->L1[0]);
		//system("pause");
	}
}

static void detslp(stec_q* stecl)
{
	FILE* fp;
	fp = fopen("F:\\CPPCode\\datas\\MW-GF-test\\Ngf-before.txt", "w+");
	char* strs;
	strs = (char*)malloc(50);
	for (int i = 0; i < eph_num; i++) { // 显示输出
		time2str(stecl->stec[i].time, strs, 4);
		fprintf(fp, "%s ", strs);
		for (int j = 0; j < NSATCMP; j++) {
			fprintf(fp, "%12.3f ", stecl->stec[i].N_gf[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	free(strs);
	strs = NULL;

	double mean_Nmw_n[NSATCMP] = { 0 }, mean_Nmw_p[NSATCMP] = { 0 };
	double mean_Ngf_n[NSATCMP] = { 0 }, mean_Ngf_p[NSATCMP] = { 0 };
	double sqr_n[NSATCMP] = { 0 }, sqr_p[NSATCMP] = { 0 };
	double ephs[NSATCMP] = { 0 }; int idx[NSATCMP] = { 0 };

	int brk = 0;
	for (int k1 = 0; k1 < eph_num; k1++) {
		for (int i = 0; i < NSATCMP; i++) {
			if (k1 >= 592 && i == 5)brk++; /* 断点测试 */
			if (stecl->stec[k1].N_mw[i] == 0 || stecl->stec[k1].N_gf[i] == 0) {
				idx[i]++; continue;
			}
			if (k1 == 0) { /* 第一个历元，初始化 */
				mean_Nmw_n[i] = mean_Nmw_p[i] = stecl->stec[k1].N_mw[i];
				mean_Ngf_n[i] = mean_Ngf_p[i] = stecl->stec[k1].N_gf[i];
				sqr_n[i] = sqr_p[i] = 0.0;
				ephs[i] = 1.0;
			}
			else { /* 如果在该历元前一个观测值都没有，或者长期失锁，重新初始化 */
				if (mean_Nmw_n[i] == 0 && mean_Nmw_p[i] == 0 && sqr_n[i] == 0 && sqr_p[i] == 0) {
					mean_Nmw_n[i] = mean_Nmw_p[i] = stecl->stec[k1].N_mw[i];
					mean_Ngf_n[i] = mean_Ngf_p[i] = stecl->stec[k1].N_gf[i];
					sqr_n[i] = sqr_p[i] = 0.0;
					ephs[i] = 1.0; idx[i] = 0;
				}
				else { /* 该历元前的观测值连续 */
					mean_Nmw_n[i] = ((ephs[i] - 1) / ephs[i]) * mean_Nmw_p[i] + (1 / ephs[i]) * stecl->stec[k1].N_mw[i];
					mean_Ngf_n[i] = ((ephs[i] - 1) / ephs[i]) * mean_Ngf_p[i] + (1 / ephs[i]) * stecl->stec[k1].N_gf[i];
					sqr_n[i] = ((ephs[i] - 1) / ephs[i]) * sqr_p[i] + (1 / ephs[i]) * pow((stecl->stec[k1].N_mw[i] - mean_Nmw_p[i]), 2);
					ephs[i]++;
				}
			}

			if (k1 > 0 && k1 < eph_num && ephs[i] > 1) { /* 当前历元非首个历元 */
				if (stecl->stec[k1 - 1].N_gf[i] != 0 && stecl->stec[k1 - 1].N_mw[i] != 0) {                       /* <1>.之前历元连续 */
					idx[i] = 0;
					double thr_mw = 5.0 * sqrt(sqr_n[i]);
					double dlt_mw = fabs(stecl->stec[k1].N_mw[i] - mean_Nmw_p[i]);
					double dlt_gf = fabs(stecl->stec[k1].N_gf[i] - stecl->stec[k1 - 1].N_gf[i]);

					if (dlt_mw > thr_mw || dlt_gf > 0.6) {                                                      /* <1.1>.连续且超过阈值 */
						//if (i == 5)printf("发现周跳—> 历元%4d  卫星C0%d\n", (k1 + 1), i + 1); /* 周跳输出 */
						if (k1 < eph_num - 1) {
							if (stecl->stec[k1 + 1].N_mw[i] != 0 && abs(stecl->stec[k1].N_mw[i] - stecl->stec[k1 + 1].N_mw[i]) >= 1) {
								if (i == 5)printf("发现连续野值—> 历元%4d  卫星C0%d\n", k1 + 1, i + 1); /* 野值输出 */
								stecl->stec[k1].L1[i] = stecl->stec[k1].L2[i] = stecl->stec[k1].P1[i] = stecl->stec[k1].P2[i] = 0;
								stecl->stec[k1].N_mw[i] = stecl->stec[k1].N_gf[i] = 0;
								ephs[i]--; idx[i]++;
								continue;
							}
							else repslp(stecl, &mean_Nmw_p, &sqr_n, k1, i, idx[i]);
						}
						else {
							repslp(stecl, &mean_Nmw_p, &sqr_n, k1, i, idx[i]);
						}
					}
				}
				else {                                                                                            /* <2>.之前历元失锁 */
					if (ephs[i] == 1) {                                                                           /* <2.1>.之前没有数据或刚被重置 */
						if (k1 < eph_num - 1) {
							if ((stecl->stec[k1 + 1].N_mw[i] != 0 && abs(stecl->stec[k1].N_mw[i] - stecl->stec[k1 + 1].N_mw[i]) >= 1) || stecl->stec[k1 + 1].N_mw[i] == 0) {
								stecl->stec[k1].L1[i] = stecl->stec[k1].L2[i] = stecl->stec[k1].P1[i] = stecl->stec[k1].P2[i] = 0;
								stecl->stec[k1].N_mw[i] = stecl->stec[k1].N_gf[i] = 0;
								mean_Nmw_n[i] = mean_Nmw_p[i] = sqr_n[i] = sqr_p[i] = 0;
								continue;
							}
						}
						else {
							stecl->stec[k1].L1[i] = stecl->stec[k1].L2[i] = stecl->stec[k1].P1[i] = stecl->stec[k1].P2[i] = 0;
							stecl->stec[k1].N_mw[i] = stecl->stec[k1].N_gf[i] = 0;
						}
						continue;
					}
					else {                                                                                        /* <2.2>.之前历元未被重置 */
						if (idx[i] > 1) {                                                                         /* <2.2.1>.失锁间隔大于2 */
							if (k1 < eph_num - 1) {
								if (stecl->stec[k1 + 1].N_mw[i] != 0 && abs(stecl->stec[k1].N_mw[i] - stecl->stec[k1 + 1].N_mw[i]) >= 1 || stecl->stec[k1 + 1].N_mw[i] == 0) {
									if (i == 5)printf("发现失锁野值—> 历元%4d  卫星C0%d\n", k1 + 1, i + 1); /* 野值输出 */
									stecl->stec[k1].L1[i] = stecl->stec[k1].L2[i] = stecl->stec[k1].P1[i] = stecl->stec[k1].P2[i] = 0;
									stecl->stec[k1].N_mw[i] = stecl->stec[k1].N_gf[i] = 0;
									mean_Nmw_n[i] = mean_Nmw_p[i] = sqr_n[i] = sqr_p[i] = 0;
									continue;
								}
								mean_Nmw_n[i] = mean_Nmw_p[i] = stecl->stec[k1].N_mw[i];
								mean_Ngf_n[i] = mean_Ngf_p[i] = stecl->stec[k1].N_gf[i];
								sqr_n[i] = sqr_p[i] = 0.0;
								ephs[i] = 1.0; idx[i] = 0;
							}
							else {
								if (i == 5)printf("发现失锁野值—> 历元%4d  卫星C0%d\n", k1 + 1, i + 1); /* 野值输出 */
								stecl->stec[k1].L1[i] = stecl->stec[k1].L2[i] = stecl->stec[k1].P1[i] = stecl->stec[k1].P2[i] = 0;
								stecl->stec[k1].N_mw[i] = stecl->stec[k1].N_gf[i] = 0;
							}
							continue;
						}
						else {                                                                                   /* <2.2.2>.之前失锁时间较短，可拼接计算 */
							double thr_mw = 5.0 * sqrt(sqr_n[i]);
							double dlt_mw = fabs(stecl->stec[k1].N_mw[i] - mean_Nmw_p[i]);
							double dlt_gf = fabs(stecl->stec[k1].N_gf[i] - stecl->stec[k1 - (int)(idx[i] + 1)].N_gf[i]);

							if (dlt_mw > thr_mw || dlt_gf > 0.6) {
								//if (i == 5)printf("发现周跳—> 历元%4d  卫星C0%d\n", k1 + 1, i + 1); /* 周跳输出 */
								if (k1 < eph_num - 1) {
									if (stecl->stec[k1 + 1].N_mw[i] != 0 && abs(stecl->stec[k1].N_mw[i] - stecl->stec[k1 + 1].N_mw[i]) >= 1) {
										if (i == 5)printf("发现失锁野值—> 历元%4d  卫星C0%d\n", k1 + 1, i + 1); /* 野值输出 */
										stecl->stec[k1].L1[i] = stecl->stec[k1].L2[i] = stecl->stec[k1].P1[i] = stecl->stec[k1].P2[i] = 0;
										stecl->stec[k1].N_mw[i] = stecl->stec[k1].N_gf[i] = 0;
										ephs[i]--; idx[i]++;
										continue;
									}
									else repslp(stecl, &mean_Nmw_p, &sqr_n, k1, i, idx[i]);
								}
								else {
									repslp(stecl, &mean_Nmw_p, &sqr_n, k1, i, idx[i]);
								}
							}
						}
					}
				}
			}
			if (k1 > 0) {
				mean_Nmw_p[i] = ((ephs[i] - 1.0) / ephs[i]) * mean_Nmw_p[i] + (1.0 / ephs[i]) * stecl->stec[k1].N_mw[i];
				mean_Ngf_p[i] = ((ephs[i] - 1.0) / ephs[i]) * mean_Ngf_p[i] + (1.0 / ephs[i]) * stecl->stec[k1].N_gf[i];
				sqr_p[i] = ((ephs[i] - 1) / ephs[i]) * sqr_p[i] + (1 / ephs[i]) * pow((stecl->stec[k1].N_mw[i] - mean_Nmw_p[i]), 2);
			}
		}
	}

	//FILE* fp;
	fp = fopen("F:\\CPPCode\\datas\\MW-GF-test\\Ngf-after.txt", "w+");
	//char* strs;
	strs = (char*)malloc(50);
	for (int i = 0; i < eph_num; i++) { // 显示输出
		time2str(stecl->stec[i].time, strs, 4);
		fprintf(fp, "%s ", strs);
		for (int j = 0; j < NSATCMP; j++) {
			fprintf(fp, "%12.3f ", stecl->stec[i].N_gf[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	free(strs);
	strs = NULL;

	system("pause");
}

static void repslp(stec_q* stecl, double* mean_Nmw_p, double* sqr_n, int k1, int i, int idx) {
	double dlt_Nmw = stecl->stec[k1].N_mw[i] - mean_Nmw_p[i];
	double dlt_Ngf = stecl->stec[k1].N_gf[i] - stecl->stec[k1 - idx - 1].N_gf[i];
	double dlt_N1 = (lam_1 * dlt_Ngf - lam_2 * dlt_Nmw) / (lam_1 - lam_2);
	double dlt_N2 = lam_1 * (dlt_Ngf - dlt_Nmw) / (lam_1 - lam_2);

	double res = 100; /* 预设残差 */
	double L1_temp, L2_temp;
	for (int j1 = -10; j1 < 11; j1++) {
		for (int j2 = -10; j2 < 11; j2++) {
			double L1 = stecl->stec[k1].L1[i];
			double L2 = stecl->stec[k1].L2[i];
			double P1 = stecl->stec[k1].P1[i];
			double P2 = stecl->stec[k1].P2[i];
			L1 -= (double)(round(dlt_N1) + j1);                          /* L1整数周修复 */
			L2 -= (double)(round(dlt_N2) + j2);                          /* L2整数周修复 */

			double L_wl = a_wl * lam_1 * L1 - b_wl * lam_2 * L2;
			double P_nl = a_nl * P1 + b_nl * P2;
			double N_mw_rep = (L_wl - P_nl) / lam_wl;            /* 修复Nmw */

			double L_gf = lam_1 * L1 - lam_2 * L2;
			double N_gf_rep = L_gf / lam_1;                      /* 修复Ngf */

			double dlt_Nmw_rep = (round(dlt_N1) + (double)j1) - (round(dlt_N2) + (double)j2);
			//double sqr_rep = ((ephs[i] - 1) / ephs[i]) * sqr_p[i] + (1 / ephs[i]) * pow((N_mw_rep - mean_Nmw_p[i]), 2);
			//if (k1 == 302 && i == 0 && j1 == 0 && j2 == 0)brk++;
			double thr_mw = 5.0 * sqrt(sqr_n[i]);
			double dlt_mw = fabs(N_mw_rep - mean_Nmw_p[i]);
			double dlt_gf = fabs(N_gf_rep - stecl->stec[k1 - idx - 1].N_gf[i]);

			if (dlt_mw > thr_mw || dlt_gf > 0.15) {        /* 再次判断修复后的周跳是否超出阈值 */
				continue;
			}

			if (fabs(dlt_Nmw_rep - dlt_Nmw) < res) {
				res = fabs(dlt_Nmw_rep - dlt_Nmw);
				L1_temp = L1;
				L2_temp = L2;
			}
		}
	}
	stecl->stec[k1].L1[i] = L1_temp;
	stecl->stec[k1].L2[i] = L2_temp;
	stecl->stec[k1].N_mw[i] = cal_MW(stecl->stec[k1].L1[i], stecl->stec[k1].L2[i], stecl->stec[k1].P1[i], stecl->stec[k1].P2[i]);
	stecl->stec[k1].N_gf[i] = cal_GF(stecl->stec[k1].L1[i], stecl->stec[k1].L2[i], stecl->stec[k1].P1[i], stecl->stec[k1].P2[i]);
}

static void calSTEC(stec_q* stecl) {
	double r = pow(FREQ1_CMP / FREQ2_CMP, 2);
	double wtk = 1.0 / (20.0 * 30.0); /* 权重系数 */
	double d_P_n, d_L_n, I_P_n, I_L_n, main_I;
	double d_P_p, d_L_p, I_P_p;
	double d_I, main_I_p;
	int ni = 0;
	double main_I_temp, d_I_temp;

	for (int x = 0; x < 2880; x++) {
		if ((double)stecl->stec[x].time.time != 0.0) ni++;
	}

	for (int i = 0; i < NSATCMP; i++) {
		int bk = 0;
		for (int j = 0, j_idx = 0; j < ni; j++) {
			if (j == 393 && i == 12)bk++;
			if (stecl->stec[j].L1[i] == 0 || stecl->stec[j].L2[i] == 0 ||
				stecl->stec[j].P1[i] == 0 || stecl->stec[j].P2[i] == 0) {
				stecl->stec[j].stec_value[i] = 0;
			}
			else {
				/* QHY -----------------------------------------------------------------------------------------------------*/
				d_P_n = stecl->stec[j].P2[i] - stecl->stec[j].P1[i]; // 双频伪距差值 P2 - P1
				I_P_n = d_P_n / (r - 1); // 当前时刻伪距电离层延迟
				d_L_n = stecl->stec[j].L1[i] * lam_1 - stecl->stec[j].L2[i] * lam_2; // 双频相位差值，需要乘以波长
				I_L_n = d_L_n / (r - 1);// 当前时刻相位电离层延迟
				main_I = 0;
				if (j == 0) {
					j_idx = j;
					main_I = I_P_n;
					stecl->stec[j].stec_value[i] = main_I / (40.3 / pow(FREQ1_CMP, 2)) / 1e16;
				}
				else {
					if (stecl->stec[j - 1].stec_value[i] != 0) {
						j_idx = j;

						d_P_p = stecl->stec[j - 1].P2[i] - stecl->stec[j - 1].P1[i];
						I_P_p = d_P_p / (r - 1); // 上一历元的L1伪距电离层延迟
						d_L_p = stecl->stec[j - 1].L1[i] * lam_1 - stecl->stec[j - 1].L2[i] * lam_2;

						main_I_p = stecl->stec[j - 1].stec_value[i] * (40.3 / pow(FREQ1_CMP, 2)) * 1e16;
						d_I = (d_L_n - d_L_p) / (r - 1); // L1相位电离层延迟变化率
						d_I_temp = d_I;
						main_I = wtk * I_P_n + (1 - wtk) * (main_I_p + d_I); // 通过权重，延迟变化率，上一时刻平滑值求解当前时刻平滑值
						main_I_temp = main_I;

						double tmp = main_I / (40.3 / pow(FREQ1_CMP, 2)) / 1e16;
						stecl->stec[j].stec_value[i] = main_I / (40.3 / pow(FREQ1_CMP, 2)) / 1e16;
					}
					else {
						if (j - j_idx <= 10) {
							main_I = wtk * I_P_n + (1 - wtk) * (main_I_temp + d_I_temp * (double)(j - j_idx));
							j_idx = j;
							stecl->stec[j].stec_value[i] = main_I / (40.3 / pow(FREQ1_CMP, 2)) / 1e16;
						}
						else {
							j_idx = j;
							main_I = I_P_n;
							stecl->stec[j].stec_value[i] = main_I / (40.3 / pow(FREQ1_CMP, 2)) / 1e16;
						}
					}
				}
				/* YJR ---------------------------------------------------------------------------------------------------*/
				//int n = 30; // 五个连续历元平滑伪距
				//if (2880 - j <= n)
				//{
				//	n = 2880 - j;
				//}
				//double sum = 0;
				//double index = 0; // 记录数据为0的个数
				//double flag = 0;

				//for (int k2 = 0; k2 < n; k2++){
				//	if (stecl->stec[j + k2].L1[i] == 0 ||stecl->stec[j + k2].L2[i] == 0 ||
				//		stecl->stec[j + k2].P1[i] == 0 ||stecl->stec[j + k2].P1[i] == 0) flag = 1;

				//	sum = sum + stecl->stec[j + k2].L1[i] * lam_1 - stecl->stec[j + k2].L2[i] * lam_2
				//		      + stecl->stec[j + k2].P1[i] - stecl->stec[j + k2].P2[i];
				//}
				//if (flag == 1)continue; // 如果flag为1，则跳过计算当前历元的平滑伪距步骤，进入下一次循环
				//stecl->stec[j].stec_value[i] = stecl->stec[j].L1[i] * lam_1 - stecl->stec[j].L2[i] * lam_2 - (1 / (n - index)) * sum; // 载波相位是周期值 要乘波长
				//stecl->stec[j].stec_value[i] = stecl->stec[j].stec_value[i] / (0.105);
				/* ---------------------------------------------------------------------------------------------------------*/
			}
		}
	}
	FILE* fp;
	fp = fopen("F:\\CPPCode\\datas\\STEC\\HD-stec.txt", "w+");
	char* strs;
	strs = (char*)malloc(50);
	for (int i = 0; i < ni; i++) { // 显示输出
		time2str(stecl->stec[i].time, strs, 4);
		fprintf(fp, "%s ", strs);
		for (int j = 0; j < NSATCMP; j++) {
			fprintf(fp, "%12.3f ", stecl->stec[i].stec_value[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	free(strs);
	strs = NULL;
	return;
}

/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          double tu        I   processing unit time (s) (0:all)
*          prcopt_t *popt   I   processing options
*          solopt_t *sopt   I   solution options
*          filopt_t *fopt   I   file options
*          char   **infile  I   input files (see below)
*          int    n         I   number of input files
*          char   *outfile  I   output file ("":stdout, see below)
*          char   *rov      I   rover id list        (separated by " ")
*          char   *base     I   base station id list (separated by " ")
* return : status (0:ok,0>:error,1:aborted)
* notes  : input files should contain observation data, navigation data, precise
*          ephemeris/clock (optional), sbas log file (optional), ssr message
*          log file (optional) and tec grid file (optional). only the first
*          observation data file in the input files is recognized as the rover
*          data.
*
*          the type of an input file is recognized by the file extension as ]
*          follows:
*              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
*              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
*              .lex,.LEX            : qzss lex message log files
*              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
*              .*i,.*I              : tec grid files (ionex)
*              .fcb,.FCB            : satellite fcb
*              others               : rinex obs, nav, gnav, hnav, qnav or clock
*
*          inputs files can include wild-cards (*). if an file includes
*          wild-cards, the wild-card expanded multiple files are used.
*
*          inputs files can include keywords. if an file includes keywords,
*          the keywords are replaced by date, time, rover id and base station
*          id and multiple session analyses run. refer reppath() for the
*          keywords.
*
*          the output file can also include keywords. if the output file does
*          not include keywords. the results of all multiple session analyses
*          are output to a single output file.
*
*          ssr corrections are valid only for forward estimation.
*-----------------------------------------------------------------------------*/
extern int postpos(gtime_t ts, gtime_t te, double ti, double tu,
	const prcopt_t* popt, const solopt_t* sopt,
	const filopt_t* fopt, char** infile, int n, char* outfile,
	const char* rov, const char* base)
{
	gtime_t tts, tte, ttte;
	double tunit, tss;
	int i, j, k, nf, stat = 0, week, flag = 1, index[MAXINFILE] = { 0 };
	char* ifile[MAXINFILE], ofile[1024], * ext;

	trace(3, "postpos : ti=%.0f tu=%.0f n=%d outfile=%s\n", ti, tu, n, outfile);

	/* open processing session */
	if (!openses(popt, sopt, fopt, &navs, &pcvss, &pcvsr)) return -1;

	if (ts.time != 0 && te.time != 0 && tu >= 0.0) {
		if (timediff(te, ts) < 0.0) { /* 判断1：如果ts>te,返回0 */
			showmsg("error : no period");
			closeses(&navs, &pcvss, &pcvsr);
			return 0;
		}
		for (i = 0; i < MAXINFILE; i++) {
			if (!(ifile[i] = (char*)malloc(1024))) {  /* 判断2：如果分配内存失败，返回-1 */
				for (; i >= 0; i--) free(ifile[i]);
				closeses(&navs, &pcvss, &pcvsr);
				return -1;
			}
		}
		if (tu == 0.0 || tu > 86400.0 * MAXPRCDAYS) tu = 86400.0 * MAXPRCDAYS;
		settspan(ts, te);
		tunit = tu < 86400.0 ? tu : 86400.0;
		tss = tunit * (int)floor(time2gpst(ts, &week) / tunit);

		for (i = 0;; i++) { /* for each periods */
			tts = gpst2time(week, tss + i * tu);
			tte = timeadd(tts, tu - DTTOL);
			if (timediff(tts, te) > 0.0) break;
			if (timediff(tts, ts) < 0.0) tts = ts;
			if (timediff(tte, te) > 0.0) tte = te;

			strcpy(proc_rov, "");
			strcpy(proc_base, "");
			if (checkbrk("reading    : %s", time_str(tts, 0))) {
				stat = 1;
				break;
			}
			for (j = k = nf = 0; j < n; j++) {
				ext = strrchr(infile[j], '.');

				if (ext && (!strcmp(ext, ".rtcm3") || !strcmp(ext, ".RTCM3"))) { /* 读RTCM文件 */
					strcpy(ifile[nf++], infile[j]);
				}
				else {
					/* include next day precise ephemeris or rinex brdc nav */
					ttte = tte;
					if (ext && (!strcmp(ext, ".sp3") || !strcmp(ext, ".SP3") ||
						!strcmp(ext, ".eph") || !strcmp(ext, ".EPH"))) {
						ttte = timeadd(ttte, 3600.0);
					}
					else if (strstr(infile[j], "brdc")) {
						ttte = timeadd(ttte, 7200.0);
					}
					nf += reppaths(infile[j], ifile + nf, MAXINFILE - nf, tts, ttte, "", "");
				}
				while (k < nf) index[k++] = j;

				if (nf >= MAXINFILE) {
					trace(2, "too many input files. trancated\n");
					break;
				}
			}
			if (!reppath(outfile, ofile, tts, "", "") && i > 0) flag = 0;

			/* execute processing session */
			stat = execses_b(tts, tte, ti, popt, sopt, fopt, flag, ifile, index, nf, ofile, rov, base);

			if (stat == 1) break;
		}
		for (i = 0; i < n && i < MAXINFILE; i++) free(ifile[i]);
	}
	else if (ts.time != 0) {
		for (i = 0; i < n && i < MAXINFILE; i++) {
			if (!(ifile[i] = (char*)malloc(1024))) {
				for (; i >= 0; i--) free(ifile[i]);
				return -1;
			}
			reppath(infile[i], ifile[i], ts, "", "");
			index[i] = i;
		}
		reppath(outfile, ofile, ts, "", "");

		/* execute processing session */
		stat = execses_b(ts, te, ti, popt, sopt, fopt, 1, ifile, index, n, ofile, rov,
			base);

		for (i = 0; i < n && i < MAXINFILE; i++) free(ifile[i]);
	}
	else {
		for (i = 0; i < n; i++) index[i] = i;

		/* execute processing session */
		stat = execses_b(ts, te, ti, popt, sopt, fopt, 1, infile, index, n, outfile, rov,
			base);
	}
	/* close processing session */
	closeses(&navs, &pcvss, &pcvsr);

	//system("pause");

	return stat;
}