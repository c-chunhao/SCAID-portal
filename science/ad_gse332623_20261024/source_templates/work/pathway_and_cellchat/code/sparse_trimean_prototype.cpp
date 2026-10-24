#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <cmath>
using namespace Rcpp;

// Finite nonnegative input only. Retain original cell indices, sort positive
// values once per gene, and account for all implicit zeros in every group.
// [[Rcpp::export]]
List scaid_prepare_positive_order(NumericMatrix transposed) {
  const int n=transposed.nrow(), genes=transposed.ncol();
  List values(genes),cells(genes);
  for(int j=0;j<genes;++j) {
    std::vector<std::pair<double,int>> positive;
    for(int i=0;i<n;++i) {
      const double x=transposed(i,j);
      if(!R_FINITE(x) || x<0)stop("Prototype accepts only verified finite nonnegative input");
      if(x>0)positive.emplace_back(x,i);
    }
    std::sort(positive.begin(),positive.end(),[](const std::pair<double,int>& a,const std::pair<double,int>& b){
      return a.first<b.first || (a.first==b.first && a.second<b.second);
    });
    NumericVector v(positive.size());IntegerVector ix(positive.size());
    for(size_t k=0;k<positive.size();++k){v[k]=positive[k].first;ix[k]=positive[k].second;}
    values[j]=v;cells[j]=ix;
  }
  return List::create(_["values"]=values,_["cells"]=cells,_["n_cells"]=n,
                      _["gene_names"]=colnames(transposed));
}

// Return the four official quartiles, not an algebraically simplified mean.
// R's installed mean() combines these in the wrapper to preserve its exact
// accumulation and correction behavior.
// [[Rcpp::export]]
NumericVector scaid_positive_group_quartiles(List prepared,IntegerVector group,int groups) {
  const int n=as<int>(prepared["n_cells"]);
  if(group.size()!=n || groups<1)stop("Invalid group axis");
  std::vector<int> size(groups,0);
  for(int i=0;i<n;++i){int g=group[i]-1;if(g<0||g>=groups)stop("Missing or invalid group code");++size[g];}
  for(int g=0;g<groups;++g)if(size[g]<1)stop("Empty group level");
  List values=prepared["values"],cells=prepared["cells"];
  const int genes=values.size();
  NumericVector result(4*genes*groups);
  result.attr("dim")=IntegerVector::create(4,genes,groups);
  for(int j=0;j<genes;++j) {
    NumericVector v=values[j];IntegerVector ix=cells[j];
    if(v.size()==0)continue;
    std::vector<int> count(groups,0),seen(groups,0),zeros(groups,0);
    for(R_xlen_t k=0;k<ix.size();++k)++count[group[ix[k]]-1];
    std::vector<int> lower(3*groups,-1),upper(3*groups,-1);
    std::vector<double> lo(3*groups,0),hi(3*groups,0),fraction(3*groups,0);
    int remaining=0;
    for(int g=0;g<groups;++g) {
      zeros[g]=size[g]-count[g];
      for(int q=0;q<3;++q) {
        const int a=q+3*g;
        const long long numerator=static_cast<long long>(size[g]-1)*(q+1);
        const int rank=static_cast<int>(numerator/4);
        const int remainder=static_cast<int>(numerator%4);
        fraction[a]=remainder/4.0;
        if(rank>=zeros[g]){lower[a]=rank-zeros[g];++remaining;}
        // If the interpolation coefficient is zero, the upper value is not
        // used; avoiding that unnecessary rank also handles a singleton.
        if(remainder!=0 && rank+1>=zeros[g]){upper[a]=rank+1-zeros[g];++remaining;}
      }
    }
    if(remaining>0)for(R_xlen_t k=0;k<v.size() && remaining>0;++k) {
      const int g=group[ix[k]]-1,position=seen[g]++;
      for(int q=0;q<3;++q) {
        const int a=q+3*g;
        if(lower[a]==position){lo[a]=v[k];--remaining;}
        if(upper[a]==position){hi[a]=v[k];--remaining;}
      }
    }
    if(remaining!=0)stop("Unresolved positive rank");
    for(int g=0;g<groups;++g) {
      double quantile[3];
      for(int q=0;q<3;++q) {
        const int a=q+3*g;
        if(fraction[a]==0)quantile[q]=lo[a];
        else {
          // Installed collapse::fquantile(type=7): a + gamma * (b-a).
          // Keep both elementary operations rounded as double, no FMA.
          volatile double difference=hi[a]-lo[a];
          volatile double product=fraction[a]*difference;
          quantile[q]=lo[a]+product;
        }
      }
      const R_xlen_t offset=4*(j+genes*g);
      result[offset]=quantile[0];result[offset+1]=quantile[1];
      result[offset+2]=quantile[1];result[offset+3]=quantile[2];
    }
  }
  return result;
}
