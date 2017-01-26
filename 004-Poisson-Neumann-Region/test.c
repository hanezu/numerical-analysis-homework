//
// Created by 朱華春 on 2016/12/23.
//
#include <stdio.h>
int N;
#define h (N+1)

int main(){
    printf("%d\n", h);
    N = 55;
    printf("%d\n", h);
    N = 45;
    printf("%d\n", h);

}