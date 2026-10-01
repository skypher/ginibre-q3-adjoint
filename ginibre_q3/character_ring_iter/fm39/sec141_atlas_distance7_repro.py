from pathlib import Path
import argparse
import sys

base = Path("ginibre_q3/character_ring_iter/fm39/sec138_coverage_atlas_repro.py")
source = base.read_text()

def once(old, new):
    global source
    count = source.count(old)
    if count != 1:
        raise SystemExit(f"patch anchor count {count}: {old[:70]!r}")
    source = source.replace(old, new, 1)

once(
    "R3_ONE_HAT_H1_SUFFIX,NR};",
    "R3_ONE_HAT_H1_SUFFIX,FM68,FM69,FM64D5,FM64ZERO,FM67,"
    "FM52B2,FM52MIN2,ADV1M53,FM70,FM71,FM77C,FM77I,NR};"
)
once(
    '"r3_hs_ht_h1a","r3_one_hatS_h1a"};',
    '"r3_hs_ht_h1a","r3_one_hatS_h1a",'
    '"FM-MECH68_distance_le_4",'
    '"FM-MECH69_multi_large_small_background",'
    '"FM-MECH64_two_core_distance_le_5",'
    '"FM-MECH64_two_core_min_ae_zero",'
    '"FM-MECH67_two_core_core_sum_ge_background",'
    '"FM-MECH52_one_core_b_le_2",'
    '"FM-MECH52_one_core_min_ae_le_1",'
    '"ADV1_M53_one_core_distance_3_5",'
    '"FM-MECH70_E_two_core_min_ae_le_6",'
    '"FM-MECH71_FLAGGED_distance_le_6",'
    '"FM-MECH77_cubic_threshold_FLAGGED",'
    '"FM-MECH77_criterion1_FLAGGED"};'
)
once(
    "hit[NR],grp[NH][41][ND],bysum[61];int minsum",
    "hit[NR],grp[NH][41][ND],bysum[61],marginal[NR],"
    "structure[2][5];int minsum"
)

predicates = r''' if(d>=0&&d<=4){if(!mask)S->marginal[FM68]++;bit(mask,FM68);}
 int c1=0,c2m=0,c2p=0,j3=0,j4=0,large=0;long long S69=0,Q69=0;
 for(int i=0;i<len;i++){int n=a[i].n,sg=a[i].s;if(n==1)c1++;if(n==2&&sg<0)c2m++;if(n==2&&sg>0)c2p++;if(n==3)j3++;if(n==4)j4++;if(n>=5){large=1;S69+=(long long)(n+1)*(n+1);if(n+1>Q69)Q69=n+1;}}
 long long A2=c1+2LL*c2m+2LL*(j3+j4),B2=A2+2LL*c2p;
 if(large&&S69>0&&((B2>=768LL*S69&&B2>=262144LL)||B2>=8192LL*S69)){if(!mask)S->marginal[FM69]++;bit(mask,FM69);}
 if(high==2&&d>=0&&d<=5){if(!mask)S->marginal[FM64D5]++;bit(mask,FM64D5);}
 if(high==2&&(suffix==0||e==0)){if(!mask)S->marginal[FM64ZERO]++;bit(mask,FM64ZERO);}
 int coreSum=0;for(int i=0;i<len;i++)if(a[i].n>=3)coreSum+=a[i].n;
 if(high==2&&coreSum>=suffix+e+2*pos2){if(!mask)S->marginal[FM67]++;bit(mask,FM67);}
 if(high==1&&pos2<=2){if(!mask)S->marginal[FM52B2]++;bit(mask,FM52B2);}
 if(high==1&&(suffix<=1||e<=1)){if(!mask)S->marginal[FM52MIN2]++;bit(mask,FM52MIN2);}
 if(high==1&&d>=3&&d<=5){if(!mask)S->marginal[ADV1M53]++;bit(mask,ADV1M53);}
 if(nh+ns==2&&pos2==0&&(suffix<=6||e<=6)){if(!mask)S->marginal[FM70]++;bit(mask,FM70);}
 if(d>=0&&d<=6){if(!mask)S->marginal[FM71]++;bit(mask,FM71);}
 if(mx>=3&&d>=0){
  int skipped=0,t77=0,b77=0,k77=0,sig377=0;
  for(int i=0;i<len;i++){int n=a[i].n,sg=a[i].s;if(!skipped&&n==mx){skipped=1;continue;}if(n==1)t77++;else if(n==2&&sg<0)t77+=2;else if(n==2&&sg>0)b77++;else if(n>=3){k77++;if(n==3)sig377+=sg;}}
  long long N77=t77+2LL*b77+k77,s77=d+1,C77=3200LL*s77*s77*s77+760LL*s77*s77;
  if(N77>=C77){if(!mask)S->marginal[FM77C]++;bit(mask,FM77C);}
  long long q77=3*N77-1140LL*s77*s77;
  if(N77>=380LL*s77*s77&&N77*q77*q77>=28800LL*sig377*sig377*s77*s77*s77){if(!mask)S->marginal[FM77I]++;bit(mask,FM77I);}
 }
'''
once(
    "for(int k=0;k<NR;k++)if(mask>>k&1)S->hit[k]++;*dp=d;return mask;",
    predicates + "for(int k=0;k<NR;k++)if(mask>>k&1)S->hit[k]++;*dp=d;return mask;"
)

classes = r'''S->uncovered++;if(d==7||d==8){int bp=0,c1x=0,c2mx=0,j3x=0,j4x=0;long long Sx=0,Qx=0;for(int z=0;z<st->len;z++){int n=st->a[z].n,sg=st->a[z].s;if(n==1)c1x++;if(n==2&&sg<0)c2mx++;if(n==2&&sg>0)bp++;if(n==3)j3x++;if(n==4)j4x++;if(n>=5){Sx+=(long long)(n+1)*(n+1);if(n+1>Qx)Qx=n+1;}}long long A2x=c1x+2LL*c2mx+2LL*(j3x+j4x),B2x=A2x+2LL*bp,Wx=Sx+16LL*j3x+25LL*j4x;int residual8=Sx>0&&Wx<(1LL<<21)*Qx*Qx&&(B2x<768LL*Sx||(A2x<512&&B2x<262144));int cl=(h==1)?0:((h==2&&bp>0)?1:(residual8?2:(mx<=5?3:4)));S->structure[d-7][cl]++;}I val=checked?chk:exact(st);'''
once("S->uncovered++;I val=checked?chk:exact(st);", classes)

once(
    "for(int i=0;i<=60;i++)z.bysum[i]+=s->bysum[i];for(int j=0;j<s->nex;j++)",
    "for(int i=0;i<=60;i++)z.bysum[i]+=s->bysum[i];"
    "for(int k=0;k<NR;k++)z.marginal[k]+=s->marginal[k];"
    "for(int x=0;x<2;x++)for(int y=0;y<5;y++)z.structure[x][y]+=s->structure[x][y];"
    "for(int j=0;j<s->nex;j++)"
)
once(
    "for(int h=0;h<NH;h++)for(int m=1;m<=40;m++){int any=0;",
    "for(int h=0;phase<0&&h<NH;h++)for(int m=1;m<=40;m++){int any=0;"
)
once(
    "for(int j=0;j<z.nex;j++){char b[512];",
    'for(int k=FM68;k<=FM77I;k++)if(z.marginal[k])printf("MARGINAL phase=%d theorem=%s count=%llu\\\\n",phase,names[k],z.marginal[k]);'
    'const char*cn[5]={"H1_one_large","two_cores_b_ge_1","MECH69_residual8","labels_le_5","other"};'
    'for(int x=0;x<2;x++)for(int y=0;y<5;y++)if(z.structure[x][y])printf("STRUCT phase=%d distance=%d class=%s count=%llu\\\\n",phase,x+7,cn[y],z.structure[x][y]);'
    "for(int j=0;j<z.nex;j++){char b[512];"
)
once("S->lists%2000000==0", "S->lists%500000==0")

ap = argparse.ArgumentParser()
ap.add_argument("--phase", choices=("0", "1", "both"), default="both")
ap.add_argument("--threads", type=int, default=1)
args = ap.parse_args()
sys.argv = [str(base), "--phase", args.phase, "--threads", str(args.threads)]
exec(compile(source, str(base), "exec"),
     {"__name__": "__main__", "__file__": str(base)})