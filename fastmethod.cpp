#include "branch_vec.hh"
/** 
 * @file fastmethod.cpp
 * @author Carl Boettiger, <cboettig@gmail.com>
 * @section DESCRIPTION
 * Fast method keeps track of the competition matrix explicitly
 * this avoids recalculating competition coefficients every time 
 * a birth or death event occurs.  This matrix is updated every
 * time a mutation occurs or a population type goes extinct.  Note  
 * that the coefficients are still the multiplied by the population 
 * sizes of each type after each birth/death event.  
 * 
 *
 * */


void print_C(vector<CRow> &cmatrix)
{
	vector<CRow>::iterator v;
	vector<double>::iterator w;
	for(v=cmatrix.begin(); v != cmatrix.end(); v++){
		for(w=v->begin(); w != v->end(); w++){
			printf("  %lf  ", *w);
		}
		printf("\n");
	}
	printf("\n");
}

void grow_C(vector<pop> &poplist, vector<CRow> &cmatrix, par_list * pars)
{
	vector<pop>::iterator p;
// vector<CRow>::iterator v = cmatrix.begin();
	int v = 0;
	double x = ( --poplist.end() )->trait;
	CRow row;
	double c;
	for(p=poplist.begin(); p != --poplist.end(); p++){
		c = C( x, p->trait, pars );
		row.push_back(c);					// add to final row
		cmatrix[v].push_back(c);			// add to last column
		v++;
	}
	row.push_back(1);						// competes against self
	cmatrix.push_back(row);					// append the final row to matrix

}


void shrink_C(vector<CRow> &cmatrix, int id)
{
	vector<CRow>::iterator v;
	vector<double>::iterator w;
	for(v =cmatrix.begin(); v != cmatrix.end(); v++){
		w = v->begin()+id;
		v->erase(w);
	}
	v = cmatrix.begin()+id;
	cmatrix.erase(v);

}


void fast_update_rates(vector<pop> &poplist, vector<CRow> &cmatrix, int id, char type, par_list * pars)
{
	vector<pop>::iterator p;
	vector<CRow>::iterator v;
	vector<double> row;
	
	vector<double>::iterator r;

	if(id >= 0){
		row = cmatrix[id];
	} else {
		v = (--cmatrix.end() ); //mutant case
		row = *v;
	}
		int j = 0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		switch(type)
		{
			case 'b' : 	
				p->compete += row[j];
				break;

			case 'd' : 	
				p->compete -= row[j];
				break;
		}
		j++;
		p->birth = R* p->popsize;
		p->death = R* p->popsize * (p->compete) * Kinv( p->trait, pars );
	}

}

double init_compete(vector<pop> &poplist, double y, par_list * pars)
{
	double compete = 0;
	vector<pop>::iterator p;
	for(p=poplist.begin(); p != poplist.end(); p++){
		compete += p->popsize * C(y, p->trait, pars);
	}
	return compete;
}


void event_and_rates(gsl_rng * rng, vector<pop> &poplist, double sum, vector<CRow> &cmatrix, par_list * pars)
{
	vector<pop>::iterator p;
	double threshhold = sum*gsl_rng_uniform(rng);
	double cumulative = 0;
	int id = 0; 
	char type='b';
	int who=0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		cumulative += p->birth;
		if( cumulative  > threshhold ){ 
			if(gsl_rng_uniform(rng) > pars->mu){
				++p->popsize;
				type = 'b';
				who = id;

			} else {
				double y = p->trait + gsl_ran_gaussian_ziggurat(rng,pars->sigma_mu);
				double c = init_compete(poplist, y, pars);
				pop pop_i = {1, y, 1., 1., p->trait, c};
				poplist.push_back(pop_i);

				grow_C(poplist, cmatrix, pars);
				type = 'b';
				who = -1;

			}
			break;
		}
		cumulative += p->death;
		if( cumulative  > threshhold ){ 
			--p->popsize;
			type = 'd';
			who = id;
			if( p->popsize == 0 ){
				poplist.erase( p );

				shrink_C(cmatrix, id); // destroy row and column of cmatrix
			}
			break;
		}
		id++;
	}
	fast_update_rates(poplist, cmatrix, who, type, pars);
}


