
#include <Rcpp.h>
using namespace Rcpp;

//' overlapsListvsList
//'
//' This function intersects a list of string vectors with another list of string vectors.
//' Outputs the number of items in the overlap between a pair of string vectors.
//' Rows are vectors from lt1, columns are vectors from lt2
//' adpated from: https://jokergoo.github.io/2023/04/05/speed-up-over-representation-enrichment-analysis/
//'
//' @param lt1 A list of vectors of strings
//' @param lt2 A list of vectors of strings
//' @return out a matrix[i,j] of integers overlaps between lt1[[i]] and lt2[[j]]
//' @export
// [[Rcpp::export]]
NumericMatrix overlapsListvsList(List lt1, List lt2) {

    int n1 = lt1.size(); // Rows
    int n2 = lt2.size(); // Columns
    NumericMatrix out( n1 , n2 );

    for(int row = 0; row < n1; row++) {
	//this gene set acts as reference
        StringVector row_v = as<StringVector>(lt1[row]);
	std::unordered_set<String> set1;
	set1.insert(row_v.begin(), row_v.end());
        for(int col = 0; col < n2; col++) {
            StringVector col_v = as<StringVector>(lt2[col]);
            int olap = 0;
            //for(int gene = 0; gene < col_v.size(); gene++) {
            for(int gene = 0; gene < col_v.size(); gene++) {
                if (set1.find(col_v[gene]) != set1.end()) {
                     olap++;
                }
            }
            out(row,col) = olap;
        }
    }
    return out;
}
