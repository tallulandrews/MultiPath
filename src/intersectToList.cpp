
#include <Rcpp.h>
using namespace Rcpp;

//' intersectToList
//'
//' This function intersects a list of string vectors with a single string vector.
//'
//' @param lt A list of vectors of strings
//' @param x A vector of strings
//' @export
// [[Rcpp::export]]
List intersectToList(List lt, StringVector x) {

    int n = lt.size();
    List out(n);

    std::unordered_set<String> seen;
    seen.insert(x.begin(), x.end());

    for(int i = 0; i < n; i++) {

        StringVector v = as<StringVector>(lt[i]);
        LogicalVector l(v.size());

        std::unordered_set<String> seen2;

        for(int j = 0; j < v.size(); j ++) {
            l[j] = seen.find(v[j]) != seen.end() && seen2.insert(v[j]).second;
        }

        out[i] = v[l];
    }

    return out;
}
