//
// Created by 朱華春 on 2016/12/23.
//
#include <stdio.h>
#include <stdlib.h>

#define N (10)
typedef struct stu{
    int name;
}stu;


int main(){
    stu **school = malloc(N * sizeof(stu));
    stu *student;
    for(int i = 0; i < N; i ++){
        student= malloc(sizeof(stu));
        student->name = i;

        school[i] = student;
        printf("%d ",school[i]->name);

    }
    printf("\n");

    for(int i = 0; i < N; i ++) {
        printf("%d ", school[i]->name);
    }
    free(student);
    for(int i = 0; i < N; i ++){
        free(school[i]);
    }
    free(school);



}