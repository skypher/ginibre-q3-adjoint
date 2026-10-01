import argparse
import ctypes
import glob
import os
import subprocess
import sys

def count_lists(limit, max_label, high_cap):
    # Coefficients of prod_n (1-z^n)^(-2), with sign parity
    # and a cap on occurrences of labels >= 3.
    dp = [[[0, 0] for _ in range(high_cap + 1)]
          for __ in range(limit + 1)]
    dp[0][0][0] = 1
    for n in range(1, max_label + 1):
        for sign in (-1, 1):
            for total in range(n, limit + 1):
                for high in range(high_cap + 1):
                    h2 = high + (n >= 3)
                    if h2 > high_cap:
                        continue
                    for parity in (0, 1):
                        dp[total][h2][parity ^ (sign < 0)] += \
                            dp[total - n][high][parity]
    return sum(dp[s][h][0]
               for s in range(1, limit + 1)
               for h in range(high_cap + 1))

assert count_lists(30, 30, 10) == 1114614
assert count_lists(60, 40, 3) == 78037215
print("Independent list counts: 1114614, 78037215")

SRC = r'''
#include <pthread.h>
#include <atomic>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
typedef __int128 I; typedef unsigned long long U;
struct Tok{int n,s;}; const int ML=60,MTH=32,NH=20,ND=31;
enum{PARITY,FUSION,CHARACTER,M41,M47,M48,M485,M38,M39,M35,M39D3,M4452,M52MIN,M52CUT,M51,M50,ADV1D3,M53D45,M53D6,OL,TWO,THREE,LL,DSAS,T3R,G0B,G0E,CENTERED,LLM,LLMS,EVENH,SQUARE,ESTRIP,EQ2,ME55,ESHORT,LS,LR4,ADV2,R1_HONLY,R1_ONE_HAT,R1_BIG_PART,R1_TWO_HAT_H3,R1_THREE_HAT_H2,R2_HONLY_T1,R2_HONLY_D1,R2_ONE_HAT4,R3_HONLY5,R3_TWO_H1_SUFFIX,R3_ONE_HAT_H1_SUFFIX,NR};
const char* names[NR]={"PARITY_ZERO","FUSION_SUPPORT_ZERO","CHARACTER_ALL_PLUS","FM-MECH41_labels_le_2","FM-MECH47_labels_le_3","FM-MECH48_labels_le_4","FM-MECH48_labels_le_5_Hge71","FM-MECH38_HACq_d_le_2","FM-MECH39_HACq_other_labels_ge_d","FM-MECH35_HAC_d3_tle3_or_vle1","FM-MECH39_HACq_d3_tle1","FM-MECH44_45_one_extra_le6","FM-MECH52_one_extra_min_ae_le1","FM-MECH52_one_extra_quadratic_cutoff","FM-MECH51_B2","FM-MECH50_Q2plus","ADV1_d3_one_extra","FM-MECH53_d4_d5_one_extra","FM-MECH53_d6_subregion","Theorem_OL","Theorem_TWO","Theorem_THREE","Theorem_LL","Theorem_DS_AS","Theorem_T3R","Theorem_G0B","Theorem_G0E","Theorem_centered_equal_label","Theorem_LLm","Theorem_LLm_S","H_only_even_parts_no_suffix","H_only_equal_square","E_strips","E_q2_plus_or_minus","FM-MECH55_E_regions","E_short_psi_arc","Theorem_LS","Theorem_LR4","ADV2_T_ge_2pow21","r1_H_only_all","r1_one_hatS_all","r1_LemmaBP_bigpart","r1_two_hatS_le3_H","r1_three_hatS_two_H","r2_H_only_le6_or_7_min4","r2_H_only_distance1","r2_one_hatS_le4_H","r3_H_only_le5","r3_hs_ht_h1a","r3_one_hatS_h1a"};
struct Ex{int len,sum,d,h,mx;Tok a[ML];I val;};
struct Stats{U lists,covered,uncovered,zero,negative,checks,hit[NR],grp[NH][41][ND],bysum[61];int minsum,nex;Ex ex[20];I minval;Ex minex;};
struct State{Tok a[ML];int len,sum,minus,high;};
struct Worker{int id,nt,phase,limit,cap,emit;};
struct Task{State st;int start,root;};
static Task tasks[100000]; static int ntasks=0; static std::atomic<int> nexttask(0);
static Stats stats[MTH];static Tok toks[80];static int ntok,emitall;
static pthread_mutex_t outlock=PTHREAD_MUTEX_INITIALIZER;
static thread_local I crow[61];static thread_local char ibuf[8][64];static thread_local int ibpos=0;
static char* istr(I x){ibpos=(ibpos+1)%8;char*p=ibuf[ibpos];int neg=x<0,k=0;if(neg)x=-x;do{p[k++]='0'+(int)(x%10);x/=10;}while(x&&k<60);if(neg)p[k++]='-';p[k]=0;for(int i=0;i<k/2;i++){char z=p[i];p[i]=p[k-1-i];p[k-1-i]=z;}return p;}
static void utc(char*b){time_t t=time(0);struct tm q;gmtime_r(&t,&q);snprintf(b,32,"%04d-%02d-%02dT%02d:%02d:%02dZ",q.tm_year+1900,q.tm_mon+1,q.tm_mday,q.tm_hour,q.tm_min,q.tm_sec);}
static I C(int n,int k){if(k<0||k>n)return 0;if(k>n-k)k=n-k;I z=1;for(int j=1;j<=k;j++)z=z*(n-k+j)/j;return z;}
static I Cat(int m){if(m<0||(m&1))return 0;int k=m/2;return C(2*k,k)/(k+1);}
static void row(int a,int e){int N=a+e;for(int k=0;k<=N;k++){I z=0;int lo=k-a;if(lo<0)lo=0;int hi=k<e?k:e;for(int j=lo;j<=hi;j++){I t=C(e,j)*C(a,k-j);z+=(j&1)?-t:t;}crow[k]=z;}}
static I cg(int j,int N){return j<0||j>N?0:crow[j];}
static I Bv(int j,int N){return cg(j-1,N)+cg(j+1,N);}
static I Wv(int j,int i,int N){return Bv(i,N)*cg(j,N)-cg(i,N)*Bv(j,N);}
static int shortarc(int j,int i,int N){for(int k=j+1;k<=i;k++)if(Wv(j,k,N)<=0)return 0;return i>j;}
static int me55(int a,int e,int j,int i,int N){int s=i-j;if(s<4)return 0;long long sig=N+2,del=a-e,X=2LL*j-N;if(4LL*(del<0?-del:del)>sig)return 0;int m=s/2;return(sig<=4*X&&2*X<=sig)||(X>=0&&2*X<=sig&&3LL*m*(X+2LL*m-2)>=16*sig);}
static int twobasic(int a,int e,int p,int q){if(p<q){int z=p;p=q;q=z;}if(q<=1)return 1;if(a<=2||e<=2||(a-e<=1&&a-e>=-1))return 1;int N=a+e,nn=N+p-q;if(nn&1)return 0;int j=nn/2,i=(N+p+q)/2+1;if(q==2)return 1;if(me55(a,e,j,i,N))return 1;row(a,e);return shortarc(j,i,N);}
static int lam(const int*v,int mask,int m){int s=0,mx=0;for(int i=0;i<m;i++)if(mask>>i&1){s+=v[i];if(v[i]>mx)mx=v[i];}int lo=2*mx-s;if(lo<0)lo=0;if((lo&1)!=(s&1))lo++;return lo;}
static int spread(const int*v,int m,int N){if(m<=1)return 1;int lim=1<<m;for(int a=1;a<lim-1;a++)if(lam(v,a,m)+lam(v,(lim-1)^a,m)<=N)return 0;return 1;}
static void bit(U&mask,int k){mask|=1ULL<<k;}
static int g0e3(const int*H,int r,int a){if(r<2)return 0;int v[3]={H[0],H[1],H[2]};for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)if(v[i]<v[j]){int z=v[i];v[i]=v[j];v[j]=z;}int A=v[0],B=v[1],C0=v[2],e=2*r-3,N=a+e,nn=N+A-B-C0;if(nn&1)return 0;int al=nn/2,be=al+C0,ga=al+B,Y=A-B+C0,m=a<e?a:e;if(al<0||al>=be||be>N)return 0;
 if(ga>=N+1){if(!(C0&1)&&m>=2&&!(m&1)&&3LL*m*Y*Y<=N-m+4)return 1;if((C0&1)&&1LL*(e+1)*(Y+1)*(Y+1)<=2LL*(a+2))return 1;if(!(C0&1)){for(int d=A-B;d<=A+B;d+=2)if(d>0&&!twobasic(a,e,d,C0))return 0;return 1;}}
 if(ga==N){if((C0&1)&&!(r&1)&&1LL*(e+1)*(Y+1)*(Y+1)<=2LL*(a+2))return 1;if(!(C0&1)&&m>=2&&!(m&1)&&!((al+m/2)&1)&&3LL*m*Y*Y<=N-m+4)return 1;}return 0;}
static I exact(const State*st){I a[61][61],b[61][61];memset(a,0,sizeof(a));a[0][0]=1;int mx=0;for(int z=0;z<st->len;z++){int n=st->a[z].n,s=st->a[z].s;memset(b,0,sizeof(b));for(int x=0;x<=mx;x++)for(int y=0;y<=mx;y++){I v=a[x][y];if(!v)continue;int lo=x-n;if(lo<0)lo=n-x;for(int k=lo;k<=x+n;k+=2)b[k][y]+=v;lo=y-n;if(lo<0)lo=n-y;for(int k=lo;k<=y+n;k+=2)b[x][k]+=s*v;}mx+=n;memcpy(a,b,sizeof(a));}return a[0][0];}
static I direct(const State*st){int total=st->sum;I a[15][15],b[15][15];memset(a,0,sizeof(a));a[0][0]=1;for(int z=0;z<st->len;z++){int n=st->a[z].n,s=st->a[z].s;I u0[16]={0},u1[16]={0},v[16]={0};u0[0]=1;if(n>=1)u1[1]=1;for(int k=2;k<=n;k++){memset(v,0,sizeof(v));for(int j=0;j<k;j++)v[j+1]+=u1[j];for(int j=0;j<=k-2;j++)v[j]-=u0[j];memcpy(u0,u1,sizeof(u0));memcpy(u1,v,sizeof(u1));}memset(b,0,sizeof(b));for(int x=0;x<=total;x++)for(int y=0;y<=total;y++)if(a[x][y])for(int j=0;j<=n;j++)if(u1[j]){if(x+j<=total)b[x+j][y]+=a[x][y]*u1[j];if(y+j<=total)b[x][y+j]+=a[x][y]*u1[j]*s;}memcpy(a,b,sizeof(a));}I ans=0;for(int x=0;x<=total;x+=2)for(int y=0;y<=total-x;y+=2)ans+=a[x][y]*Cat(x)*Cat(y);return ans;}
static U coverage(const State*st,int mx,int high,Stats*S,int*dp){const Tok*a=st->a;int len=st->len,sum=st->sum,minus=st->minus;U mask=0;int dn=sum-2*mx,d=-999;
 if(sum&1)bit(mask,PARITY);else if(dn<0)bit(mask,FUSION);else{d=dn/2;if(d<=2)bit(mask,M38);}
 if(d>=0){int rem=0,ok=1;for(int i=0;i<len;i++){if(!rem&&a[i].n==mx){rem=1;continue;}if(a[i].n<d){ok=0;break;}}if(ok)bit(mask,M39);int t=0,v=0;for(int i=0;i<len;i++){t+=a[i].n==1;v+=a[i].n==2;}if(d==3){if(t<=3||v<=1)bit(mask,M35);if(t<=1)bit(mask,M39D3);}}
 int allplus=1,le2=1,le3=1,le4=1,le5=1,pos1=0,neg2=0,pos2=0,H[61],SS[61],nh=0,ns=0,highn=0,weights=0;
 for(int i=0;i<len;i++){int n=a[i].n,s=a[i].s;if(s<0)allplus=0;if(n>2)le2=0;if(n>3)le3=0;if(n>4)le4=0;if(n>5)le5=0;pos1+=(n==1&&s>0);neg2+=(n==2&&s<0);pos2+=(n==2&&s>0);if(n>=3&&s<0)H[nh++]=n;if(n>=2&&s>0)SS[ns++]=n;if(n>=3){highn++;weights+=(n+1)*(n+1);}}
 int r=minus/2,suffix=pos1+neg2,e=2*r-nh;if(allplus)bit(mask,CHARACTER);if(le2)bit(mask,M41);if(le3)bit(mask,M47);if(le4)bit(mask,M48);if(le5&&highn>=71)bit(mask,M485);
 if(high==1){int p=0,negp=0;for(int i=0;i<len;i++)if(a[i].n>=3){p=a[i].n;negp=a[i].s<0;}int ee=2*r-(negp?1:0),d1=dn>=0&&!(dn&1)?dn/2:-999;if(p<=6)bit(mask,M4452);int eps=negp?1:0,K=(p&1)?4*p*p-2:(negp?2*p*p-2:2*p*p+2);if(suffix<=1||ee<=1)bit(mask,M52MIN);if(2*r+suffix+eps+2*pos2>=2*K)bit(mask,M52CUT);if(pos2==2)bit(mask,M51);if(pos2==1)bit(mask,M50);if(pos2==0)bit(mask,OL);if(d1==3)bit(mask,ADV1D3);if(d1==4||d1==5)bit(mask,M53D45);if(d1==6&&suffix>=2&&ee>=2&&pos2>=6&&2LL*(suffix+ee)>=11LL*pos2-1)bit(mask,M53D6);}
 int hfac=nh+suffix,minpart=61,sumhat=0;for(int i=0;i<nh;i++)if(H[i]-1<minpart)minpart=H[i]-1;for(int i=0;i<ns;i++)sumhat+=SS[i];if(suffix>0)minpart=1;
 if(r==1&&ns==0)bit(mask,R1_HONLY);if(r==1&&ns==1)bit(mask,R1_ONE_HAT);
 if(r==1){int big=0;for(int i=0;i<nh;i++)if(H[i]-1>=sumhat-1)big=1;if(suffix>0&&1>=sumhat-1)big=1;if(big)bit(mask,R1_BIG_PART);}
 if(r==1&&ns==2&&hfac<=3)bit(mask,R1_TWO_HAT_H3);if(r==1&&ns==3&&hfac==2)bit(mask,R1_THREE_HAT_H2);
 if(r==2&&ns==0&&(hfac<=6||(hfac==7&&minpart<=4)))bit(mask,R2_HONLY_T1);if(r==2&&ns==0&&d==1)bit(mask,R2_HONLY_D1);if(r==2&&ns==1&&hfac<=4)bit(mask,R2_ONE_HAT4);
 if(r==3&&ns==0&&hfac<=5)bit(mask,R3_HONLY5);if(r==3&&ns==0&&nh==2)bit(mask,R3_TWO_H1_SUFFIX);if(r==3&&ns==1&&nh==0)bit(mask,R3_ONE_HAT_H1_SUFFIX);
 if(nh+ns<=1)bit(mask,OL);if(nh==2&&ns==0&&suffix==0&&r>=1)bit(mask,TWO);if(nh==3&&ns==0&&suffix==0)bit(mask,THREE);
 if(nh==3&&ns==0){int u[3]={H[0]-1,H[1]-1,H[2]-1},v[3]={H[0],H[1],H[2]};for(int i=0;i<3;i++)for(int j=i+1;j<3;j++){if(u[i]<u[j]){int z=u[i];u[i]=u[j];u[j]=z;}if(v[i]<v[j]){int z=v[i];v[i]=v[j];v[j]=z;}}
 if(u[2]>=suffix+2*r-4||u[0]-u[1]+u[2]>=suffix+2*r-3)bit(mask,LL);if(r>=2&&(suffix==2*r-3||suffix==2*r-4||suffix==2*r-2))bit(mask,DSAS);
 if(r>=2){int odds=(u[0]&1)+(u[1]&1)+(u[2]&1);if(suffix==2*r-3&&odds==1)bit(mask,T3R);if(suffix==2*r-3&&odds==3&&u[1]%4==1&&u[2]%4==1)bit(mask,T3R);
 if(u[0]+u[1]-u[2]>=suffix+2*r-3&&((suffix>=2*r-4&&suffix<=2*r-2)||(suffix<2*r-3?suffix:2*r-3)<=2))bit(mask,T3R);
 if((suffix-2*r+2)*(suffix-2*r+2)<=8*r-9){int A=v[0],B=v[1],C0=v[2],N=suffix+2*r-3,nn=N+A-B-C0;if(!(nn&1)&&nn/2+B>=N+1)bit(mask,G0B);}
 if(g0e3(H,r,suffix))bit(mask,G0E);if(u[0]==u[1]&&2LL*u[0]>=suffix+2LL*r+u[2]-2)bit(mask,CENTERED);}}
 if(nh>0&&ns==0&&suffix==0){int ok=1;for(int i=0;i<nh;i++)if((H[i]-1)&1)ok=0;if(ok)bit(mask,EVENH);}
 if(nh>=2&&ns==0&&!(nh&1)){int ok=1;for(int i=1;i<nh;i++)if(H[i]!=H[0])ok=0;if(ok)bit(mask,SQUARE);}
 int core[61],m=0;for(int i=0;i<nh;i++)core[m++]=H[i];for(int i=0;i<ns;i++)core[m++]=SS[i];int N=suffix+2*r-nh;if(nh<=2*r&&spread(H,nh,N))bit(mask,LLM);if(m<=2*r&&spread(core,m,N))bit(mask,LLMS);
 long long L15=0;for(int i=0;i<nh;i++){long long k=H[i]-1;L15+=3*k*(k+4);}for(int i=0;i<ns;i++){long long k=SS[i];L15+=5*k*(k+2);}if(15LL*(suffix+2LL*r+4)>=(2LL*r+3)*L15)bit(mask,LS);
 long long L6=0;int tw=0;for(int i=0;i<nh;i++){long long k=H[i]-1;tw+=k&1;if(!(k&1))L6+=4*k*k*(k+1)*(k+3);else L6+=(k-1)*(k-1)*(k+2)*(k+3);}for(int i=0;i<ns;i++){long long k=SS[i];tw+=k&1;if(!(k&1))L6+=12*k*k;else{long long z=k-1;L6+=4*z*z*k*(k+2);}}long long lrN=2LL*r+suffix+tw;if(!((suffix+tw)&1)&&6*(lrN+6)>=7*L6)bit(mask,LR4);
 if(nh+ns==2){int p0[2],q0=0;for(int i=0;i<nh;i++)p0[q0++]=H[i];for(int i=0;i<ns;i++)p0[q0++]=SS[i];int ee=2*r-nh;if(ee>=0){int pp=p0[0],qq=p0[1];if(pp<qq){int z=pp;pp=qq;qq=z;}int NN=suffix+ee,num=NN+pp-qq;if(!(num&1)){int j=num/2,ii=(NN+pp+qq)/2+1;if(suffix<=2||ee<=2||(suffix-ee<=1&&suffix-ee>=-1))bit(mask,ESTRIP);if(qq==2)bit(mask,EQ2);if(me55(suffix,ee,j,ii,NN))bit(mask,ME55);row(suffix,ee);if(shortarc(j,ii,NN))bit(mask,ESHORT);}}}
 int K=mx+1;if(mx>=3&&(long long)weights>=(1LL<<21)*K*K)bit(mask,ADV2);
 for(int k=0;k<NR;k++)if(mask>>k&1)S->hit[k]++;*dp=d;return mask;}
static void word(const Tok*a,int n,char*b,int cap){int k=snprintf(b,cap,"(");for(int i=0;i<n&&k<cap-20;i++)k+=snprintf(b+k,cap-k,"%s%c%d",i?",":"",a[i].s<0?'-':'+',a[i].n);snprintf(b+k,cap-k,")");}
static void addex(Stats*S,const State*st,int d,int h,int mx,I val){if(st->sum<S->minsum){S->minsum=st->sum;S->nex=0;}if(st->sum==S->minsum&&S->nex<20){Ex&e=S->ex[S->nex++];e.len=st->len;e.sum=st->sum;e.d=d;e.h=h;e.mx=mx;e.val=val;memcpy(e.a,st->a,st->len*sizeof(Tok));}if(val>0&&val<S->minval){S->minval=val;Ex&e=S->minex;e.len=st->len;e.sum=st->sum;e.d=d;e.h=h;e.mx=mx;e.val=val;memcpy(e.a,st->a,st->len*sizeof(Tok));}}
static void candidate(Stats*S,State*st,int phase,int emit,int tid){S->lists++;int mx=0,h=0;for(int i=0;i<st->len;i++){if(st->a[i].n>mx)mx=st->a[i].n;h+=st->a[i].n>=3;}
 I chk=0;int checked=0;if(S->checks<16&&st->sum<=12){chk=exact(st);I other=direct(st);if(chk!=other){printf("CHECK_FAIL exact=%s direct=%s\n",istr(chk),istr(other));abort();}S->checks++;checked=1;}
 int d=-999;U mask=coverage(st,mx,h,S,&d);if(mask){S->covered++;return;}S->uncovered++;I val=checked?chk:exact(st);if(!val)S->zero++;if(val<0)S->negative++;if(d>=0&&d<ND&&h<NH&&mx<=40)S->grp[h][mx][d]++;if(st->sum<=60)S->bysum[st->sum]++;addex(S,st,d,h,mx,val);
 if(emit){char b[512];word(st->a,st->len,b,sizeof(b));pthread_mutex_lock(&outlock);printf("UNCOVERED phase=%d sum=%d high=%d max=%d distance=%d value=%s word=%s\n",phase,st->sum,h,mx,d,istr(val),b);fflush(stdout);pthread_mutex_unlock(&outlock);}}
static void walk(Stats*S,State*st,int phase,int limit,int cap,int start,int root,int tid){if(st->len&&!(st->minus&1)){candidate(S,st,phase,emitall,tid);if(S->lists%2000000==0){char t[32];utc(t);printf("PROGRESS %s phase=%d thread=%d candidates=%llu depth=%d partial_sum=%d root=%d\n",t,phase,tid,S->lists,st->len,st->sum,root);fflush(stdout);}}
 for(int z=start;z<ntok;z++){Tok q=toks[z];if(st->sum+q.n>limit)continue;int hh=st->high+(q.n>=3);if(cap>=0&&hh>cap)continue;st->a[st->len++]=q;st->sum+=q.n;st->minus+=(q.s<0);st->high=hh;walk(S,st,phase,limit,cap,z,root,tid);st->high-=(q.n>=3);st->minus-=(q.s<0);st->sum-=q.n;st->len--;}}
static void* worker(void*arg){Worker*w=(Worker*)arg;Stats*S=&stats[w->id];for(;;){int z=nexttask.fetch_add(1);if(z>=ntasks)break;Task*q=&tasks[z];walk(S,&q->st,w->phase,w->limit,w->cap,q->start,q->root,w->id);}return 0;}
static void report(int phase,int threads){Stats z;memset(&z,0,sizeof(z));z.minsum=1000;z.minval=((I)1<<126);for(int t=0;t<threads;t++){Stats*s=&stats[t];z.lists+=s->lists;z.covered+=s->covered;z.uncovered+=s->uncovered;z.zero+=s->zero;z.negative+=s->negative;z.checks+=s->checks;for(int k=0;k<NR;k++)z.hit[k]+=s->hit[k];for(int h=0;h<NH;h++)for(int m=0;m<=40;m++)for(int d=0;d<ND;d++)z.grp[h][m][d]+=s->grp[h][m][d];for(int i=0;i<=60;i++)z.bysum[i]+=s->bysum[i];for(int j=0;j<s->nex;j++){Ex e=s->ex[j];if(e.sum<z.minsum){z.minsum=e.sum;z.nex=0;}if(e.sum==z.minsum&&z.nex<20)z.ex[z.nex++]=e;}if(s->minval<z.minval){z.minval=s->minval;z.minex=s->minex;}}
 char t[32];utc(t);printf("RESULT %s phase=%d lists=%llu covered=%llu uncovered=%llu value_zero=%llu value_negative=%llu crosschecks=%llu min_uncovered_sum=%d min_positive_value=%s\n",t,phase,z.lists,z.covered,z.uncovered,z.zero,z.negative,z.checks,z.minsum,istr(z.minval));for(int k=0;k<NR;k++)if(z.hit[k])printf("HIT phase=%d theorem=%s count=%llu\n",phase,names[k],z.hit[k]);for(int h=0;h<NH;h++)for(int m=1;m<=40;m++){int any=0;for(int d=0;d<ND;d++)any|=(z.grp[h][m][d]!=0);if(any){printf("GROUP phase=%d high=%d max=%d",phase,h,m);for(int d=0;d<ND;d++)if(z.grp[h][m][d])printf(" %d:%llu",d,z.grp[h][m][d]);printf("\n");}}for(int j=0;j<z.nex;j++){char b[512];word(z.ex[j].a,z.ex[j].len,b,sizeof(b));printf("SMALLEST phase=%d rank=%d sum=%d high=%d max=%d distance=%d value=%s word=%s\n",phase,j+1,z.ex[j].sum,z.ex[j].h,z.ex[j].mx,z.ex[j].d,istr(z.ex[j].val),b);}if(z.minval<((I)1<<126)){char b[512];word(z.minex.a,z.minex.len,b,sizeof(b));printf("MINVALUE phase=%d value=%s sum=%d word=%s\n",phase,istr(z.minex.val),z.minex.sum,b);}printf("UNCOVERED_BY_SUM phase=%d",phase);for(int i=0;i<=60;i++)if(z.bysum[i])printf(" %d:%llu",i,z.bysum[i]);printf("\n");fflush(stdout);}
extern "C" const char* fm138_theorem_name(int id){return id>=0&&id<NR?names[id]:"";}
extern "C" U fm138_membership(int len,const int*labels,const int*signs){if(len<0||len>ML)return 0;State st={};st.len=len;for(int i=0;i<len;i++){int n=labels[i],s=signs[i];if(n<1||(s!=1&&s!=-1)||st.sum+n>60)return 0;st.a[i]={n,s};st.sum+=n;st.minus+=(s<0);st.high+=(n>=3);}if(st.minus&1)return 0;int mx=0;for(int i=0;i<len;i++)if(st.a[i].n>mx)mx=st.a[i].n;if(!((st.sum<=30&&mx<=30)||(st.sum<=60&&mx<=40&&st.high<=3)))return 0;Stats dummy={};int d;return coverage(&st,mx,st.high,&dummy,&d);}
extern "C" void run_atlas(int phase,int threads,int emit){if(threads<1)threads=1;if(threads>MTH)threads=MTH;emitall=emit;int maxn=phase?40:30,limit=phase?60:30,cap=phase?3:-1;ntok=0;for(int n=1;n<=maxn;n++){toks[ntok++]={n,-1};toks[ntok++]={n,1};}for(int j=0;j<threads;j++){memset(&stats[j],0,sizeof(stats[j]));stats[j].minsum=1000;stats[j].minval=((I)1<<126);}ntasks=0;nexttask.store(0);State st;for(int i=0;i<ntok;i++){Tok x=toks[i];if(x.n>limit)continue;st.len=1;st.a[0]=x;st.sum=x.n;st.minus=x.s<0;st.high=x.n>=3;if(cap<0||st.high<=cap){if(!(st.minus&1))candidate(&stats[0],&st,phase,emit,0);for(int j=i;j<ntok;j++){Tok y=toks[j];if(st.sum+y.n>limit)continue;int hh=st.high+(y.n>=3);if(cap>=0&&hh>cap)continue;State pair=st;pair.a[1]=y;pair.len=2;pair.sum+=y.n;pair.minus+=(y.s<0);pair.high=hh;if(!(pair.minus&1))candidate(&stats[0],&pair,phase,emit,0);for(int k=j;k<ntok;k++){Tok z=toks[k];if(pair.sum+z.n>limit)continue;int h3=pair.high+(z.n>=3);if(cap>=0&&h3>cap)continue;if(ntasks>=100000)abort();Task&q=tasks[ntasks++];q.st=pair;q.st.a[2]=z;q.st.len=3;q.st.sum+=z.n;q.st.minus+=(z.s<0);q.st.high=h3;q.start=k;q.root=i;}}}}char t[32];utc(t);printf("START %s phase=%d labels<=%d total<=%d high_cap=%d threads=%d tasks=%d\n",t,phase,maxn,limit,cap,threads,ntasks);fflush(stdout);pthread_t th[MTH];Worker ws[MTH];for(int j=0;j<threads;j++){ws[j]={j,threads,phase,limit,cap,emit};pthread_create(&th[j],0,worker,&ws[j]);}for(int j=0;j<threads;j++)pthread_join(th[j],0);report(phase,threads);}
'''

def memfd(name, data=b""):
    fd = os.memfd_create(name, 0)
    if data:
        os.write(fd, data)
    return fd

ap = argparse.ArgumentParser(description="Exact FM-SEC138 coverage atlas")
ap.add_argument("--phase", choices=["0", "1", "both"], default="both")
ap.add_argument("--threads", type=int, default=32)
ap.add_argument("--emit-all", action="store_true",
                help="print each uncovered signed list and exact integer value")
ap.add_argument("--match",
                help="one signed list as comma-separated tokens, e.g. -1,-1,+1,+1")
args = ap.parse_args()

# Independent exact cardinality recurrence.
print("Independent list counts: 1114614, 78037215")

sfd = memfd("fm138.cpp", SRC.encode())
cc = subprocess.run(
    ["g++", "-std=c++17", "-O3", "-fPIC", "-pthread", "-S",
     "-x", "c++", "-o", "-", f"/proc/self/fd/{sfd}"],
    pass_fds=(sfd,), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if cc.returncode:
    sys.stderr.buffer.write(cc.stderr)
    sys.exit(cc.returncode)
afd = memfd("fm138.s", cc.stdout)
ofd = memfd("fm138.o")
asr = subprocess.run(
    ["as", "--64", "-o", f"/proc/self/fd/{ofd}", f"/proc/self/fd/{afd}"],
    pass_fds=(afd, ofd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if asr.returncode:
    sys.stderr.buffer.write(asr.stderr)
    sys.exit(asr.returncode)
sofd = memfd("fm138.so")
gccs = sorted(glob.glob("/usr/lib/gcc/x86_64-linux-gnu/*/libgcc.a"))
if not gccs:
    raise RuntimeError("libgcc.a was not found")
ld = subprocess.run(
    ["ld.gold", "-shared", "-o", f"/proc/self/fd/{sofd}",
     f"/proc/self/fd/{ofd}", gccs[-1], "-lpthread"],
    pass_fds=(ofd, sofd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if ld.returncode:
    sys.stderr.buffer.write(ld.stderr)
    sys.exit(ld.returncode)

lib = ctypes.CDLL(f"/proc/self/fd/{sofd}")
if args.match is not None:
    parts = args.match.split(",")
    labels, signs = [], []
    for token in parts:
        if len(token) < 2 or token[0] not in "+-" or not token[1:].isdigit():
            raise SystemExit("tokens must look like +3 or -4")
        signs.append(-1 if token[0] == "-" else 1)
        labels.append(int(token[1:]))
    if sum(s < 0 for s in signs) % 2:
        raise SystemExit("the FM3 input needs an even number of minus signs")
    aa = (ctypes.c_int * len(labels))(*labels)
    ss = (ctypes.c_int * len(signs))(*signs)
    lib.fm138_membership.argtypes = (
        ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
    lib.fm138_membership.restype = ctypes.c_ulonglong
    lib.fm138_theorem_name.argtypes = (ctypes.c_int,)
    lib.fm138_theorem_name.restype = ctypes.c_char_p
    mask = lib.fm138_membership(len(labels), aa, ss)
    print([lib.fm138_theorem_name(i).decode()
           for i in range(64) if mask >> i & 1])
else:
    run = lib.run_atlas
    run.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_int)
    run.restype = None
    phases = (0, 1) if args.phase == "both" else (int(args.phase),)
    for phase in phases:
        run(phase, max(1, min(args.threads, 32)), int(args.emit_all))

