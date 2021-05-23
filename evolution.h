#ifndef __evolution_H__
#define __evolution_H__


#define DNA_size 23			//23个基因 
#define pop_size 1000		//100个数 

extern double mutation_rate;			//变异率 
extern int generation;					//进行几代 
extern double pop_decimal[pop_size];				//十进制基因种群 
extern int pop_binary[pop_size][DNA_size];		//转化为二进制基因种群 
extern double fitness_val[pop_size];			//适应度值 
extern double pop_val[pop_size+1];				//个体的概率值 

extern void first_generation(double*);
extern void get_fitness(double*,double*,double*);
extern void translate_to_binary(double*,int (*)[DNA_size]);
extern void translate_to_decimal(double*,int (*)[DNA_size]);  
extern void crossover(int (*)[DNA_size]);
extern int select_individual(double*); 
extern int select_DNA();
extern void mutate(int (*)[DNA_size]);
extern int find_maxval(double*);
extern void print_maxval(); 
double convert_num(double);

#endif
