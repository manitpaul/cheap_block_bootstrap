#include <Rcpp.h>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

// Four fixed prognostic strata, two response categories, two treatment arms.
double contrast(const double *c0,const double *c1) {
  double N0=0,N1=0,U=0;
  for(int j=0;j<8;++j){N0+=c0[j];N1+=c1[j];}
  if(N0<=0||N1<=0)return NA_REAL;
  for(int h=0;h<4;++h){double t=c1[2*h]+c1[2*h+1];
    double n=t+c0[2*h]+c0[2*h+1];
    if(n>0)U+=c1[2*h+1]-t*(c0[2*h+1]+c1[2*h+1])/n;
  }
  return (N0+N1)*U/(N0*N1);
}
void multi(const double *pr,int number,double *out) {
  double sum=0;for(int j=0;j<8;++j)sum+=pr[j];
  int left=number;double mass=sum;
  for(int j=0;j<7;++j){double p=mass>0?pr[j]/mass:0;
    p=std::max(0.0,std::min(1.0,p));
    out[j]=left>0?R::rbinom(left,p):0;left-=int(out[j]);mass-=pr[j];}
  out[7]=left;
}
void hyper(const double *pr,int number,double *out) {
  double mass=0;for(int j=0;j<8;++j)mass+=pr[j];int left=number;
  for(int j=0;j<7;++j){out[j]=left>0?R::rhyper(pr[j],mass-pr[j],left):0;
    left-=int(out[j]);mass-=pr[j];}
  out[7]=left;
}
std::vector<int> expand(const double *pr) {
  std::vector<int>x;x.reserve(int(std::accumulate(pr,pr+8,0.0)));
  for(int j=0;j<8;++j)for(int k=0;k<int(pr[j]);++k)x.push_back(j);return x;
}
void shuffle_rng(std::vector<int>& x){for(int i=int(x.size())-1;i>0;--i){int j=int(R::runif(0,i+1));std::swap(x[i],x[j]);}}
// The generic BLB kernel draws multinomial weights over b actual support points.
void weighted_generic(const std::vector<int>& bag,int n,double *out){
  std::fill(out,out+8,0.0);int rem=n;
  for(int j=0;j<int(bag.size())-1;++j){int w=R::rbinom(rem,1.0/(bag.size()-j));out[bag[j]]+=w;rem-=w;}
  out[bag.back()]+=rem;
}
NumericVector summarize(std::vector<double>& errors) {
  int B=errors.size();double mean=std::accumulate(errors.begin(),errors.end(),0.0)/B,v=0;
  for(double x:errors)v+=(x-mean)*(x-mean);
  std::sort(errors.begin(),errors.end());
  int lo=int(std::floor((B+1)*.025+1e-10));
  int hi=lo+int(std::ceil((B+1)*.95-1e-10));
  if(lo<1||hi>B)stop("B too small for requested basic ranks");
  return NumericVector::create(_["se"]=std::sqrt(v/(B-1)),_["qlo"]=errors[lo-1],_["qhi"]=errors[hi-1]);
}

// [[Rcpp::export]]
double ist_estimate(NumericMatrix counts) {
  double c0[8],c1[8];for(int j=0;j<8;++j){c0[j]=counts(0,j);c1[j]=counts(1,j);}return contrast(c0,c1);
}

// [[Rcpp::export]]
NumericVector ist_resample(NumericMatrix counts,int B,std::string kind,int size,bool generic=false) {
  RNGScope scope;double c0[8],c1[8],t0[8],t1[8];int n0=0,n1=0;
  for(int j=0;j<8;++j){c0[j]=counts(0,j);c1[j]=counts(1,j);n0+=int(c0[j]);n1+=int(c1[j]);}
  int n=n0+n1;double theta=contrast(c0,c1);std::vector<double>err(B);
  std::vector<int>x0,x1;if(generic){x0=expand(c0);x1=expand(c1);}
  int d0=std::max(1,std::min(n0,int(std::round(size*double(n0)/n))));
  int d1=std::max(1,std::min(n1,size-d0));
  if(kind=="ordinary"||kind=="mn") {
    if(kind=="ordinary"){d0=n0;d1=n1;}
    for(int r=0;r<B;++r){
      if(generic){std::fill(t0,t0+8,0.0);std::fill(t1,t1+8,0.0);
        for(int i=0;i<d0;++i)t0[x0[int(R::runif(0,n0))]]++;
        for(int i=0;i<d1;++i)t1[x1[int(R::runif(0,n1))]]++;
      }else{multi(c0,d0,t0);multi(c1,d1,t1);}
      for(int j=0;j<8;++j){t0[j]*=double(n0)/d0;t1[j]*=double(n1)/d1;}
      err[r]=(contrast(t0,t1)-theta)*(kind=="mn"?std::sqrt(double(d0+d1)/n):1);
    }
    NumericVector ans=summarize(err);return NumericVector::create(_["se"]=ans[0],_["qlo"]=ans[1],_["qhi"]=ans[2],_["K"]=1,_["inner_fits"]=B);
  }
  if(kind=="cbb") {
    std::vector<double>bins0(d0*8,0),bins1(d1*8,0);double remain0[8],remain1[8];
    std::copy(c0,c0+8,remain0);std::copy(c1,c1+8,remain1);
    if(generic){shuffle_rng(x0);shuffle_rng(x1);}
    int used=0;
    for(int k=0;k<d0;++k){int number=n0/d0+(k<n0%d0);
      if(generic){for(int j=0;j<number;++j)bins0[8*k+x0[used++]]++;}
      else{hyper(remain0,number,&bins0[8*k]);for(int j=0;j<8;++j)remain0[j]-=bins0[8*k+j];}}
    used=0;
    for(int k=0;k<d1;++k){int number=n1/d1+(k<n1%d1);
      if(generic){for(int j=0;j<number;++j)bins1[8*k+x1[used++]]++;}
      else{hyper(remain1,number,&bins1[8*k]);for(int j=0;j<8;++j)remain1[j]-=bins1[8*k+j];}}
    for(int r=0;r<B;++r){std::fill(t0,t0+8,0.0);std::fill(t1,t1+8,0.0);
      for(int k=0;k<d0;++k){int j=int(R::runif(0,d0));for(int v=0;v<8;++v)t0[v]+=bins0[8*j+v];}
      for(int k=0;k<d1;++k){int j=int(R::runif(0,d1));for(int v=0;v<8;++v)t1[v]+=bins1[8*j+v];}
      double mass0=std::accumulate(t0,t0+8,0.0),mass1=std::accumulate(t1,t1+8,0.0);
      for(int j=0;j<8;++j){t0[j]*=n0/mass0;t1[j]*=n1/mass1;}
      err[r]=contrast(t0,t1)-theta;
    }
    NumericVector ans=summarize(err);return NumericVector::create(_["se"]=ans[0],_["qlo"]=ans[1],_["qhi"]=ans[2],_["K"]=1,_["inner_fits"]=B);
  }
  if(kind=="blb") {
    int K=int(std::ceil(double(n)/size));double means[3]={0,0,0},bag0[8],bag1[8],base0[8],base1[8];
    if(K>std::min(n0,n1))stop("BLB bags must contain both arms");
    double remain0[8],remain1[8];std::copy(c0,c0+8,remain0);std::copy(c1,c1+8,remain1);
    int used0=0,used1=0;if(generic){shuffle_rng(x0);shuffle_rng(x1);}
    for(int k=0;k<K;++k){
      d0=n0/K+(k<n0%K);d1=n1/K+(k<n1%K);
      if(generic){std::fill(bag0,bag0+8,0);std::fill(bag1,bag1+8,0);
        for(int j=0;j<d0;++j)bag0[x0[used0++]]++;
        for(int j=0;j<d1;++j)bag1[x1[used1++]]++;
      }else{hyper(remain0,d0,bag0);hyper(remain1,d1,bag1);}
      for(int j=0;j<8;++j){remain0[j]-=bag0[j];remain1[j]-=bag1[j];}
      for(int j=0;j<8;++j){base0[j]=bag0[j]*n0/d0;base1[j]=bag1[j]*n1/d1;}
      double base=contrast(base0,base1);std::vector<int>bb0,bb1;
      if(generic){bb0=expand(bag0);bb1=expand(bag1);}
      for(int r=0;r<B;++r){
        if(generic){weighted_generic(bb0,n0,t0);weighted_generic(bb1,n1,t1);}else{multi(bag0,n0,t0);multi(bag1,n1,t1);}
        err[r]=contrast(t0,t1)-base;
      }
      NumericVector one=summarize(err);for(int j=0;j<3;++j)means[j]+=one[j]/K;
    }
    return NumericVector::create(_["se"]=means[0],_["qlo"]=means[1],_["qhi"]=means[2],_["K"]=K,_["inner_fits"]=double(K)*B);
  }
  stop("Unknown method");return NumericVector();
}

// [[Rcpp::export]]
NumericVector ist_reference(NumericMatrix counts,int B) {
  RNGScope scope;double c0[8],c1[8],u0[8],u1[8];int n0=0,n1=0;
  for(int j=0;j<8;++j){c0[j]=counts(0,j);c1[j]=counts(1,j);n0+=int(c0[j]);n1+=int(c1[j]);}
  NumericVector z(B);for(int r=0;r<B;++r){multi(c0,n0,u0);multi(c1,n1,u1);z[r]=contrast(u0,u1);}return z;
}

// Exact iid finite-n sampling distribution for a categorical population.
// Conditional arm probabilities sum to one; arm size itself is Binomial(n,.5).
// [[Rcpp::export]]
NumericVector ist_population_reference(NumericMatrix probabilities,int n,int reps) {
  RNGScope scope;double p0[8],p1[8],c0[8],c1[8];for(int j=0;j<8;++j){p0[j]=probabilities(0,j);p1[j]=probabilities(1,j);}
  NumericVector ans(reps);for(int r=0;r<reps;++r){int n1=R::rbinom(n,.5);multi(p0,n-n1,c0);multi(p1,n1,c1);ans[r]=contrast(c0,c1);}return ans;
}
