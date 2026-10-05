// Full MatrixMarket coordinate stream audit, plain or gzip. Does not densify.
#include <zlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <cerrno>
#include <cctype>
#include <cstdint>
static bool line(gzFile f,std::string& s){
 s.clear();char b[65536];
 while(gzgets(f,b,sizeof(b))){s+=b;if(!s.empty()&&s.back()=='\n'){s.pop_back();if(!s.empty()&&s.back()=='\r')s.pop_back();return true;}}
 if(!gzeof(f))throw std::runtime_error("gzip/decode error");return !s.empty();
}
static void whitespace(const char*&p){while(*p&&std::isspace(static_cast<unsigned char>(*p)))++p;}
static uint64_t integer(const char*&p){
 whitespace(p);if(!*p||*p=='-')throw std::runtime_error("Invalid nonnegative integer coordinate/dimension");
 errno=0;char* end;unsigned long long n=std::strtoull(p,&end,10);
 if(end==p||errno==ERANGE||(*end&&!std::isspace(static_cast<unsigned char>(*end))))throw std::runtime_error("Coordinate/dimension is not an exact integer token");
 p=end;return n;
}
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("Usage: audit_matrixmarket_gz FILE.mtx[.gz]");
 gzFile f=gzopen(argv[1],"rb");if(!f)throw std::runtime_error("open failed");gzbuffer(f,8388608);
 std::string s;if(!line(f,s))throw std::runtime_error("Empty MatrixMarket");
 std::string magic,type,storage,field,sym,extra;std::istringstream head(s);head>>magic>>type>>storage>>field>>sym;
 auto lower=[](std::string x){std::transform(x.begin(),x.end(),x.begin(),[](unsigned char c){return std::tolower(c);});return x;};
 magic=lower(magic);type=lower(type);storage=lower(storage);field=lower(field);sym=lower(sym);
 if(magic!="%%matrixmarket"||type!="matrix"||storage!="coordinate"||(field!="integer"&&field!="real")||sym!="general"||(head>>extra))throw std::runtime_error("Require MatrixMarket coordinate integer/real general header");
 do{if(!line(f,s))throw std::runtime_error("Missing dimensions");}while(s.empty()||s[0]=='%');
 const char* q=s.c_str();uint64_t nr=integer(q),nc=integer(q),declared=integer(q);whitespace(q);if(*q||nr==0||nc==0)throw std::runtime_error("Invalid dimensions");
 uint64_t n=0,nonzero=0,nonint=0,negative=0,nonfinite=0,zeros=0;double maxv=0;uint64_t first_bad_entry=0;double first_bad_value=0;
 while(line(f,s)){
  const char*p=s.c_str();whitespace(p);if(!*p||*p=='%')continue;
  uint64_t i=integer(p),j=integer(p);if(i<1||i>nr||j<1||j>nc)throw std::runtime_error("Coordinate out of declared bounds");
  whitespace(p);if(!*p)throw std::runtime_error("Missing numeric value");
  errno=0;char* end;double v=std::strtod(p,&end);if(end==p)throw std::runtime_error("Malformed numeric value");
  p=end;whitespace(p);if(*p)throw std::runtime_error("Extra/malformed MatrixMarket tokens");
  ++n;if(v!=0)++nonzero;else ++zeros;
  if(!std::isfinite(v))++nonfinite;
  else {if(v<0)++negative;if(v!=std::floor(v)){++nonint;if(first_bad_entry==0){first_bad_entry=n;first_bad_value=v;}}maxv=std::max(maxv,v);}
 }
 gzclose(f);if(n!=declared)throw std::runtime_error("Observed entry count differs from declared nnz");
 std::cout<<std::setprecision(17)<<"{\n\"format\": \"MatrixMarket\",\n\"schema_pass\": true,\n\"features\": "<<nr<<",\n\"cells\": "<<nc<<",\n\"nnz\": "<<n<<",\n\"declared_nnz\": "<<declared<<",\n\"numeric_entries_scanned\": "<<n<<",\n\"nonzero_entries\": "<<nonzero<<",\n\"explicit_zero_entries\": "<<zeros<<",\n\"noninteger_values\": "<<nonint<<",\n\"negative_values\": "<<negative<<",\n\"nonfinite_values\": "<<nonfinite<<",\n\"maximum_value\": "<<maxv<<",\n\"first_noninteger_entry_1based\": "<<first_bad_entry<<",\n\"first_noninteger_value\": "<<first_bad_value<<",\n\"full_numeric\": \""<<(nonint==0&&negative==0&&nonfinite==0?"integer_nonnegative":"noninteger_or_invalid")<<"\",\n\"matrix_market_field\": \""<<field<<"\",\n\"coordinate_integer_and_bounds_pass\": true,\n\"bounded_reader\": \"gzip/plain streaming; 8MiB compressed buffer and one coordinate line; no densification\",\n\"duplicate_coordinate_validation\": \"not performed; repeated entries may be valid MatrixMarket sums and require canonicalization before analysis\"\n}\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
