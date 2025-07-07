#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//符号のパラーメータの指定。通常[N,K,T]として、
//Nは符号の長さ、Kが符号の次元、Tは訂正エラー数
//を表す。ここではDは符号長にしている。
#define N 2039 // 符号長
#define M N //有限体の元の数
#define K (16) //符号の次元
#define E (11)    //拡大体のビット数
#define DEG N //(K * E)
#define T (K / 2) //エラーの数


unsigned short mat[N][N] = {0};
unsigned short ma[N][N] = {0};
