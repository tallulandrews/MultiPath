
#include <Rcpp.h>
using namespace Rcpp;

//' intersectToList
//'
//' This function intersects a list of string vectors with a single string vector.
//' adpated from: https://jokergoo.github.io/2023/04/05/speed-up-over-representation-enrichment-analysis/
//' originally created by: Zuguang Gu
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

//' charSort
//'
//' This function sorts a character vector
//' adpated from: https://stackoverflow.com/questions/38654019/rcpp-sort-descending
//'
//' @param x A vector of strings
//' @return The vector sorted alphabetically a->z
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector charSort(Rcpp::CharacterVector x) {
    Rcpp::CharacterVector res = Rcpp::clone(x);
    res.sort(false);
    return res;
}
