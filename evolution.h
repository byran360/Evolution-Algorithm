#ifndef __evolution_H__
#define __evolution_H__


#define DNA_size 23			//23������ 
#define pop_size 1000		//100���� 

extern double mutation_rate;			//������ 
extern int generation;					//���м��� 
extern double pop_decimal[pop_size];				//ʮ���ƻ�����Ⱥ 
extern int pop_binary[pop_size][DNA_size];		//ת��Ϊ�����ƻ�����Ⱥ 
extern double fitness_val[pop_size];			//��Ӧ��ֵ 
extern double pop_val[pop_size+1];				//����ĸ���ֵ 

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
