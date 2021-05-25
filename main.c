// This project uses evolutionary algorithms to solve the question of maximum-solving.
// If it solves the question of minimum-solving, only need I change the fitness function.

#include <stdio.h>
#include <stdlib.h>
#include "evolution.h"
	
#define DNA_size 23			//23个基因 
#define pop_size 1000		//1000个数 

double mutation_rate = 0.003;			//变异率 
int generation = 5000;					//进行几代 
double pop_decimal[pop_size];			//十进制基因种群 
int pop_binary[pop_size][DNA_size];		//转化为二进制基因种群 
double fitness_val[pop_size];			//适应度值 
double pop_val[pop_size+1];				//个体的概率值 


void main(int argc, char *argv[]) {
	srand((unsigned)time(NULL));
	int i;
	first_generation(pop_decimal);					//产生十进制的种群 
	translate_to_binary(pop_decimal,pop_binary);	//复制并转化一个二进制种群
	get_fitness(pop_decimal,pop_val,fitness_val);	//得到适应值
	for(i = 0; i < generation; i++){
		printf("\nIt is %d times circle:\n",i+1);
		crossover(pop_binary); 
		mutate(pop_binary);
		translate_to_decimal(pop_decimal,pop_binary);
		get_fitness(pop_decimal,pop_val,fitness_val);	//得到适应值
		translate_to_binary(pop_decimal,pop_binary);	//复制并转化一个二进制种群
	}
	
	print_maxval();
}
