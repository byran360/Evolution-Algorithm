

#include <stdio.h>
#include <stdlib.h>
#include "evolution.h"
	
#define DNA_size 23			//23������ 
#define pop_size 1000		//1000���� 

double mutation_rate = 0.003;			//������ 
int generation = 5000;					//���м��� 
double pop_decimal[pop_size];			//ʮ���ƻ�����Ⱥ 
int pop_binary[pop_size][DNA_size];		//ת��Ϊ�����ƻ�����Ⱥ 
double fitness_val[pop_size];			//��Ӧ��ֵ 
double pop_val[pop_size+1];				//����ĸ���ֵ 


void main(int argc, char *argv[]) {
	srand((unsigned)time(NULL));
	int i;
	first_generation(pop_decimal);					//����ʮ���Ƶ���Ⱥ 
	translate_to_binary(pop_decimal,pop_binary);	//���Ʋ�ת��һ����������Ⱥ
	get_fitness(pop_decimal,pop_val,fitness_val);	//�õ���Ӧֵ
	for(i = 0; i < generation; i++){
		printf("\nIt is %d times circle:\n",i+1);
		crossover(pop_binary); 
		mutate(pop_binary);
		translate_to_decimal(pop_decimal,pop_binary);
		get_fitness(pop_decimal,pop_val,fitness_val);	//�õ���Ӧֵ
		translate_to_binary(pop_decimal,pop_binary);	//���Ʋ�ת��һ����������Ⱥ
	}
	
	print_maxval();
}
