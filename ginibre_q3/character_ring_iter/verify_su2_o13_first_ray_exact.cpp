#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_rational;
using Rat = cpp_rational;

Rat decimal(const std::string& text) {
    bool negative = false;
    std::size_t pos = 0;
    if (!text.empty() && text[0] == '-') { negative = true; pos = 1; }
    cpp_int numerator = 0;
    cpp_int denominator = 1;
    bool after_point = false;
    for (; pos < text.size(); ++pos) {
        const char c = text[pos];
        if (c == '.') { after_point = true; continue; }
        if (c < '0' || c > '9') throw std::runtime_error("bad decimal");
        numerator = numerator * 10 + static_cast<unsigned>(c - '0');
        if (after_point) denominator *= 10;
    }
    if (negative) numerator = -numerator;
    return Rat(numerator) / Rat(denominator);
}

struct Interval {
    Rat lo;
    Rat hi;
    Interval() : lo(0), hi(0) {}
    explicit Interval(const Rat& x) : lo(x), hi(x) {}
    Interval(const Rat& a, const Rat& b) : lo(a), hi(b) { if (hi < lo) throw std::runtime_error("bad interval"); }
};
Interval operator+(const Interval& a,const Interval& b){return Interval(a.lo+b.lo,a.hi+b.hi);}
Interval operator-(const Interval& a){return Interval(-a.hi,-a.lo);}
Interval operator-(const Interval& a,const Interval& b){return a+(-b);}
Interval operator*(const Interval& a,const Interval& b){std::array<Rat,4> v{a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};return Interval(*std::min_element(v.begin(),v.end()),*std::max_element(v.begin(),v.end()));}
Interval operator*(const Interval& a,const Rat& b){return a*Interval(b);}
Interval operator/(const Interval& a,const Rat& b){return a*Interval(Rat(1)/b);}
Interval power(Interval a,unsigned n){Interval r(Rat(1));while(n){if(n&1U)r=r*a;n>>=1U;if(n)a=a*a;}return r;}
Interval square(const Interval& a){if(a.lo <= 0 && a.hi >= 0) return Interval(Rat(0),std::max(a.lo*a.lo,a.hi*a.hi));return Interval(std::min(a.lo*a.lo,a.hi*a.hi),std::max(a.lo*a.lo,a.hi*a.hi));}
Interval absolute(const Interval& a){if(a.lo>0)return a;if(a.hi<0)return -a;throw std::runtime_error("interval crosses zero");}
int sign_of(const Interval& a){if(a.lo>0)return 1;if(a.hi<0)return -1;throw std::runtime_error("unknown sign");}

using Poly = std::vector<Rat>;
void trim(Poly& p){while(p.size()>1 && p.back()==0)p.pop_back();}
Poly derivative(const Poly& p){Poly d;for(std::size_t i=1;i<p.size();++i)d.push_back(p[i]*Rat(i));if(d.empty())d.push_back(0);return d;}
std::pair<Poly,Poly> divrem(Poly a,const Poly& b){trim(a);Poly q(a.size()>=b.size()?a.size()-b.size()+1:1,Rat(0));while(!(a.size()==1&&a[0]==0)&&a.size()>=b.size()){std::size_t k=a.size()-b.size();Rat c=a.back()/b.back();q[k]=c;for(std::size_t i=0;i<b.size();++i)a[i+k]-=c*b[i];trim(a);}return{q,a};}
Rat evaluate(const Poly& p,const Rat& x){Rat y=0;for(auto it=p.rbegin();it!=p.rend();++it)y=y*x+*it;return y;}
std::vector<Poly> sturm(const Poly& p){std::vector<Poly>s{p,derivative(p)};while(!(s.back().size()==1&&s.back()[0]==0)){auto r=divrem(s[s.size()-2],s.back()).second;if(r.size()==1&&r[0]==0)break;for(Rat&x:r)x=-x;s.push_back(r);}return s;}
int variations(const std::vector<Poly>& s,const Rat& x){int last=0,v=0;for(const Poly&p:s){Rat z=evaluate(p,x);int sg=z>0?1:(z<0?-1:0);if(sg==0)continue;if(last&&sg!=last)++v;last=sg;}return v;}

Interval b4(const Interval& x){return power(x,4)-power(x,3)*Rat(3)+x*Rat(3);}
struct Term { std::string pair; Interval coefficient; Interval lambda; };
const Term& find(const std::vector<Term>& terms,const char* pair){for(const Term&t:terms)if(t.pair==pair)return t;throw std::runtime_error(std::string("missing pair ")+pair);}
struct Allocation { const char* negative; std::vector<std::pair<const char*,unsigned>> positive; };
const std::vector<Allocation> allocations{
 {"12",{{"36",100}}},
 {"13",{{"02",4},{"03",1},{"04",18},{"05",3},{"34",1},{"35",71},{"45",2}}},
 {"14",{{"02",1},{"06",55},{"34",1},{"35",32},{"36",11}}},
 {"15",{{"05",1},{"06",68},{"36",31}}},
 {"16",{{"05",59},{"06",41}}},
 {"23",{{"35",100}}},
 {"24",{{"03",95},{"35",5}}},
 {"26",{{"06",100}}},
 {"46",{{"35",100}}},
 {"56",{{"01",19},{"03",8},{"05",73}}},
};

std::vector<int> outputs(int a,int b){std::vector<int>z;for(int c=std::abs(a-b);c<=std::min(a+b,13-a-b);++c)z.push_back(c);return z;}
using State=std::map<std::pair<int,int>,cpp_int>;
State step(const State&s,int label,int sign){State r;for(const auto&kv:s){int a=kv.first.first,b=kv.first.second;for(int x:outputs(a,label))r[{x,b}]+=kv.second;for(int x:outputs(b,label))r[{a,x}]+=sign*kv.second;}for(auto it=r.begin();it!=r.end();)if(it->second==0)it=r.erase(it);else ++it;return r;}
cpp_int corner(int p){State s{{{0,0},cpp_int(1)}};s=step(s,4,-1);for(int k=0;k<1+2*p;++k)s=step(s,1,-1);auto it=s.find({0,0});return it==s.end()?cpp_int(0):it->second;}
}

int main(){try{
    const Poly polynomial{Rat(0),Rat(5),Rat(0),Rat(-15),Rat(5),Rat(9),Rat(-6),Rat(1)};
    const auto sequence=sturm(polynomial);
    const std::array<std::pair<const char*,const char*>,7> bounds{{
      {"2.82709091528520179100425514397063435588162091875494","2.82709091528520179100425514397063435588162091875495"},
      {"2.33826121271771642765254666137356094719916643791959","2.33826121271771642765254666137356094719916643791960"},
      {"1.61803398874989484820458683436563811772030917980576","1.61803398874989484820458683436563811772030917980577"},
      {"0.79094307346469305720033169039500376183868826105081","0.79094307346469305720033169039500376183868826105082"},
      {"0","0"},
      {"-0.61803398874989484820458683436563811772030917980577","-0.61803398874989484820458683436563811772030917980576"},
      {"-0.95629520146761127585713349573919906491947561772536","-0.95629520146761127585713349573919906491947561772535"}
    }};
    std::array<Interval,7> roots;
    for(std::size_t i=0;i<bounds.size();++i){Rat lo=decimal(bounds[i].first),hi=decimal(bounds[i].second);roots[i]=Interval(lo,hi);if(lo==hi){if(evaluate(polynomial,lo)!=0)throw std::runtime_error("bad exact root");}else if(variations(sequence,lo)-variations(sequence,hi)!=1)throw std::runtime_error("root isolation failed");}
    if(variations(sequence,Rat(-4))-variations(sequence,Rat(4))!=7)throw std::runtime_error("root count failed");
    std::array<Interval,7> values4,weights;for(std::size_t i=0;i<7;++i){values4[i]=b4(roots[i]);weights[i]=(Interval(Rat(3))-roots[i])/Rat(15);}
    std::vector<Term> positive,negative;
    const Poly golden{Rat(-1),Rat(-1),Rat(1)};
    if(evaluate(golden,roots[2].lo)*evaluate(golden,roots[2].hi)>0 || evaluate(golden,roots[5].lo)*evaluate(golden,roots[5].hi)>0)throw std::runtime_error("golden root check failed");
    for(std::size_t i=0;i<7;++i)for(std::size_t j=i+1;j<7;++j){if(i==2&&j==5)continue;Interval d1=roots[i]-roots[j],d4=values4[i]-values4[j];int sign=0;try{sign=sign_of(d1)*sign_of(d4);}catch(...){throw std::runtime_error(std::string("unknown sign pair ")+std::to_string(i)+std::to_string(j));}Interval lambda=square(d1);Interval coefficient=weights[i]*weights[j]*Rat(2)*absolute(d1)*absolute(d4)*power(lambda,2);Term t{std::to_string(i)+std::to_string(j),coefficient,lambda};(sign>0?positive:negative).push_back(t);}
    if(positive.size()!=10||negative.size()!=10)throw std::runtime_error("term count failed");
    std::map<std::string,Interval> used;
    for(const Allocation&a:allocations){const Term&n=find(negative,a.negative);unsigned sum=0;Interval product(Rat(1));for(const auto&e:a.positive){const Term&p=find(positive,e.first);sum+=e.second;product=product*power(p.lambda,e.second);Interval contribution=n.coefficient*Rat(e.second)/Rat(100);auto it=used.find(e.first);if(it==used.end())used.emplace(e.first,contribution);else it->second=it->second+contribution;}if(sum!=100)throw std::runtime_error("weights do not sum");if(product.lo<power(n.lambda,100).hi)throw std::runtime_error(std::string("geometric inequality failed ")+a.negative);}
    for(const Term&p:positive){auto it=used.find(p.pair);Interval amount=it==used.end()?Interval(Rat(0)):it->second;if(amount.hi>p.coefficient.lo)throw std::runtime_error("capacity failed "+p.pair);}
    if(corner(0)!=0||corner(1)!=0||corner(2)!=8)throw std::runtime_error("leaf check failed");
    std::cout<<"SU2_O13_FIRST_RAY_EXACT PASS roots=7 spectral_pairs=21 positives=10 negatives=10 floor=2 denominator=100 leaves=2 leaf_values=0,0 tail_first=8\n";
    return 0;
}catch(const std::exception&e){std::cerr<<"error: "<<e.what()<<'\n';return 1;}}
