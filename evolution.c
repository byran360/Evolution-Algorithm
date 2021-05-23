#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h> 
#include <string.h>
#include "evolution.h"

#define DNA_size 23			//23个基因 
#define pop_size 1000		//100个数 
#define range_max 5.000000	//量程 
#define range_min 0.000000

extern double mutation_rate;			//变异率 
extern int generation;					//进行几代 
extern double pop_decimal[pop_size];				//十进制基因种群 
extern int pop_binary[pop_size][DNA_size];		//转化为二进制基因种群 
extern double fitness_val[pop_size];			//适应度值 
extern double pop_val[pop_size+1];				//个体的概率值 


/*****产生初代种群*******/
void first_generation(double* pop_decimal){
	double sum;
	int i, j, a;
	for(i = 0; i < pop_size; i++){
		for(j = 0, sum = 0; j < 7; j++){
			if(j == 0)    a = rand() % 5;
			else    a = rand() % 10;
			sum = sum + a / pow(10,j);
		}
		if(sum == 0)	sum = convert_num(sum);
		pop_decimal[i] = sum;
	}

	printf("pop_decimal:\n");	
	for(i = 0; i < pop_size; i++){	
		printf("pop_decimal[%d]=%f\n",i,pop_decimal[i]);	
	}

} 

/******获取特征值*******/ 
void get_fitness(double* pop_decimal,double* pop_val,double* fitness_val){
	int n;
	double sum_fit = 0; 
	for(n = 0; n < pop_size; n++){
		fitness_val[n] = exp(pop_decimal[n]) / pop_decimal[n] + pop_decimal[n] / exp(pop_decimal[n]);
		sum_fit += fitness_val[n];
	}
	for(n = 0, pop_val[0] = 0; n < pop_size; n++)
		pop_val[n+1] = pop_val[n] + fitness_val[n] / sum_fit;  		//[0,1]概率范围展示 

	printf("\nsum_fit:%lf\n",sum_fit);
	
	printf("\nfitness_val:\n");	
	for(n = 0; n < pop_size; n++){	
		printf("fitness_val[%d]=%lf\n",n,fitness_val[n]);
	}
	
	printf("\npop_val:\n");	
	for(n = 0; n <= pop_size; n++){	
		printf("pop_val[%d]=%lf\n",n,pop_val[n]);	
	}

}

/*******翻译DNA********/ 
//十进制转换成二进制 
void translate_to_binary(double* pop_decimal,int (*pop_binary)[DNA_size]){
	int i, j;
	int arr[pop_size];
	for(i = 0; i < pop_size; i++)		//复制数组，并成为整数 
		arr[i] = pop_decimal[i] * 1000000;
	for(i = 0; i < pop_size; i++)
		for (j = DNA_size-1; j >= 0; arr[i] /= 2, j--) 
			pop_binary[i][j] = arr[i] % 2;
		
	printf("\npop_binary:");
	for(i = 0; i < pop_size; i++){
		printf("\npop_binary[%d]=",i);	
		for(j = 0; j < DNA_size; j++)
			printf("%d",pop_binary[i][j]);	
	}

}

//二进制转换成十进制 
void translate_to_decimal(double* pop_decimal,int (*pop_binary)[DNA_size]){
	int i, j, temp, mul;
	int arr[pop_size] = {0};
	double div = 1000000.0;
 	for(i = 0; i < pop_size; i++){
 		for(j = 0,mul = 22; j < DNA_size; j++, mul--) 
        	arr[i] = arr[i] + pop_binary[i][j] * pow(2,mul);
        //printf("\narr[%d]=%d\n",i,arr[i]);
	}     
	for(i = 0; i < pop_size; i++){		//复制数组，并成为整数 
		pop_decimal[i] = arr[i] / div;
		if(pop_decimal[i] > range_max)    pop_decimal[i] = convert_num(pop_decimal[i]);
	}
		
	printf("\npop_decimal_c:\n");	
	for(i = 0; i < pop_size; i++)
		printf("pop_decimal_c[%d]=%lf\n",i,pop_decimal[i]);


}


//杂交 
void crossover(int (*pop_binary)[DNA_size]){
	int cross_i, cross_j;
	int cross_dad, cross_mom, cross_cir;  			//父母 
	int change_positin;				//染色体上第几位
	int temp, pos; 
	int temp_binary[pop_size][DNA_size]; 
	for(cross_i = 0; cross_i < pop_size; cross_i++)
		for(cross_j = 0; cross_j < DNA_size; cross_j++)
			temp_binary[cross_i][cross_j] = pop_binary[cross_i][cross_j];
			
	for(cross_cir = 0, pos = 0; cross_cir < 500; cross_cir++){
		cross_dad = select_individual(pop_val);
		//printf("\ncross_dad=%d\n",cross_dad);
		cross_mom = select_individual(pop_val);
		//printf("\ncross_mom=%d\n",cross_mom);
		change_positin = select_DNA();
		//printf("\nchange_positin=%d\n",change_positin);
	//交换染色体 
		temp = temp_binary[cross_dad][change_positin];
		temp_binary[cross_dad][change_positin] = temp_binary[cross_mom][change_positin];
		temp_binary[cross_mom][change_positin] = temp;
		
		for(cross_j = 0; cross_j < DNA_size; cross_j++)
			pop_binary[pos][cross_j] = temp_binary[cross_dad][cross_j];
		pos++;
		for(cross_j = 0; cross_j < DNA_size; cross_j++)
			pop_binary[pos][cross_j] = temp_binary[cross_mom][cross_j];
		pos++;
	} 
	
	int i, j; 
	printf("\npop_binary:");
	for(i = 0; i < pop_size; i++){
		printf("\npop_binary_c[%d]=",i);	
		for(j = 0; j < DNA_size; j++)
			printf("%d",pop_binary[i][j]);	
	}
	
}

//选择一个个体作为父 
int select_individual(double* pop_val){
	double a, select_p = 0;
	int i, j;
	for(j = 1; j < 7; j++){
		a = rand() % 10;
		select_p = select_p + a / pow(10,j);
	}
	//printf("\np=%lf",p);
	for(i = 0; i <= pop_size; i++)
		if(pop_val[i] <= select_p && pop_val[i+1] > select_p)    
			return i;
	printf("It seems somewhere is wrong.");
	return -1;
}

//选择一个染色体 
int select_DNA(){
	int q;
	q = rand() % 23;		 		 //获取0~22的随机数 
	return q;
}

//变异 
void mutate(int (*pop_binary)[DNA_size]){	
	double mutation_a, mutation_p = 0;
	int mutation_self, mutation_j, mutation_posotion;
	for(mutation_j = 1; mutation_j < 4; mutation_j++){
		mutation_a = rand() % 10;
		mutation_p = mutation_p + mutation_a / pow(10,mutation_j);
	} 
	if(mutation_p < mutation_rate){
		mutation_self = select_individual(pop_val);
		mutation_posotion = select_DNA();
		if(pop_binary[mutation_self][mutation_posotion] = 0)    
			pop_binary[mutation_self][mutation_posotion] = 1;
		else	pop_binary[mutation_self][mutation_posotion] = 0;
	}
}

//找到适应值最大的数的序号 
int find_maxval(double* fitness_val){
	int max_val, i;
	for(i = 1; i < pop_size; i++)
		if(fitness_val[max_val] < fitness_val[i])
			max_val = i;		
	return max_val;
}

//打印出结果 
void print_maxval(){
	int n;
	n = find_maxval(fitness_val);
	printf("\n\nwhen x = %f, value achieve maxmimum:%f",pop_decimal[n],fitness_val[n]);
}


//回到量程内 
double convert_num(double num){
	int j, a;
	double sum;
	while(num > 5 || num <= 0){
		for(j = 0, sum = 0, a = 0; j < 7; j++){
			if(j == 0)    a = rand() % 5;
			else    a = rand() % 10;
			sum = sum + a / pow(10,j);
		}
		if(num <= range_min)	num += sum;
		if(num > range_max)		num -= sum;	
	}
	return num;	
}

