from pathlib import Path
import argparse
import sys

ap = argparse.ArgumentParser(description="Read-only FM-SEC140 exact atlas verifier")
ap.add_argument("--phase", choices=["0", "1", "both"], default="both")
ap.add_argument("--threads", type=int, default=32)
ap.add_argument("--emit-all", action="store_true")
opt = ap.parse_args()

path = Path("ginibre_q3/character_ring_iter/fm39/sec138_coverage_atlas_repro.py")
text = path.read_text()

def patch(old, new):
    global text
    assert text.count(old) == 1, (old, text.count(old))
    text = text.replace(old, new, 1)

patch(
    "R3_ONE_HAT_H1_SUFFIX,NR};",
    "R3_ONE_HAT_H1_SUFFIX,FM68,FM6369,FM64D5,FM64STRIP,FM67,"
    "FM52B2,FM52MINNEW,ADV1M53D35,FM52CUTNEW,FM70,ADV2NEW,NR};"
)
patch(
    '"r3_one_hatS_h1a"};',
    '"r3_one_hatS_h1a","FM-MECH68_distance_le_4",'
    '"FM-MECH63_69_large_labels_small_background",'
    '"FM-MECH64_two_core_distance_le_5",'
    '"FM-MECH64_two_core_min_ae_zero",'
    '"FM-MECH67_two_core_core_sum_ge_background",'
    '"FM-MECH52_one_core_b_le_2",'
    '"FM-MECH52_one_core_min_ae_le_1",'
    '"ADV1_M53_one_core_distance_3_5",'
    '"FM-MECH52_one_core_uniform_cutoff",'
    '"FM-MECH70_E_two_core_min_ae_le_6",'
    '"ADV2_weighted_cutoff"};'
)
patch("bysum[61];int minsum", "bysum[61];U marginal[11];int minsum")
patch(
    "for(int k=0;k<NR;k++)z.hit[k]+=s->hit[k];",
    "for(int k=0;k<NR;k++)z.hit[k]+=s->hit[k];"
    "for(int k=0;k<11;k++)z.marginal[k]+=s->marginal[k];"
)
patch(
    "for(int h=0;h<NH;h++)for(int m=1;m<=40;m++){int any=0;",
    r'''for(int k=0;k<11;k++)if(z.marginal[k])
 printf("MARGINAL phase=%d theorem=%s count=%llu\n",
        phase,names[NR-11+k],z.marginal[k]);
 for(int h=0;h<NH;h++)for(int m=1;m<=40;m++){int any=0;'''
)
patch(
    "int K=mx+1;if(mx>=3&&(long long)weights>=(1LL<<21)*K*K)bit(mask,ADV2);",
    r'''int K=mx+1;if(mx>=3&&(long long)weights>=(1LL<<21)*K*K)bit(mask,ADV2);
 U oldmask=mask;
 int c1p=0,c1m=0,c2p=0,c2m=0,c3=0,c4=0,largeS=0;
 int coreSum=0,bgSum=0,eps=0,p=0;
 for(int z=0;z<len;z++){
   int n=a[z].n,sg=a[z].s;
   if(n==1){bgSum+=n;if(sg>0)c1p++;else c1m++;}
   else if(n==2){bgSum+=n;if(sg>0)c2p++;else c2m++;}
   else if(n>=3){
     coreSum+=n;if(n>p)p=n;if(sg<0)eps=1;
     if(n>=5)largeS+=(n+1)*(n+1);
     else if(n==3||n==4){c3+=(n==3);c4+=(n==4);}
   }
 }
 int ae=c1p+c2m,ee=2*r-nh,mine=ae<ee?ae:ee;
 if(d>=0&&d<=4)bit(mask,FM68);
 long long B2=(long long)c1p+c1m+2LL*c2m+2LL*c2p+2LL*(c3+c4);
 if(largeS&&B2>=768LL*largeS&&B2>=262144)bit(mask,FM6369);
 if(high==2&&d>=0&&d<=5)bit(mask,FM64D5);
 if(high==2&&mine==0)bit(mask,FM64STRIP);
 if(high==2&&coreSum>=bgSum)bit(mask,FM67);
 if(high==1&&c2p<=2)bit(mask,FM52B2);
 if(high==1&&mine<=1)bit(mask,FM52MINNEW);
 if(high==1&&d>=3&&d<=5)bit(mask,ADV1M53D35);
 if(high==1&&2LL*r+ae+eps+2LL*c2p>=8LL*(p+eps)*(p+eps)-4)
   bit(mask,FM52CUTNEW);
 if(high==2&&c2p==0&&mine<=6)bit(mask,FM70);
 if(mx>=3&&(long long)weights>=(1LL<<21)*K*K)bit(mask,ADV2NEW);
 int fresh[11]={FM68,FM6369,FM64D5,FM64STRIP,FM67,FM52B2,
   FM52MINNEW,ADV1M53D35,FM52CUTNEW,FM70,ADV2NEW};
 if(!oldmask)for(int z=0;z<11;z++)
   if(mask>>fresh[z]&1){S->marginal[z]++;break;}'''
)

sys.argv = [
    str(path), "--phase", opt.phase, "--threads", str(opt.threads)
]
if opt.emit_all:
    sys.argv.append("--emit-all")
exec(compile(text, str(path), "exec"),
     {"__name__": "__main__", "__file__": str(path)})
