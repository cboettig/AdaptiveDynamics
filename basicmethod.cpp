
void update_rates(vector<pop> &poplist)
{
	vector<pop>::iterator p, q;
	for(p=poplist.begin(); p != poplist.end(); p++){
		double compete = 0;
		for(q=poplist.begin(); q != poplist.end(); q++){
			compete += q->popsize * C( q->trait, p->trait );
		}

		p->birth = R* p->popsize;
		p->death = R* p->popsize * compete/K( p->trait );
	}

}

void event(gsl_rng * rng, vector<pop> &poplist, double sum)
{
	vector<pop>::iterator p;
	double threshhold = sum*gsl_rng_uniform(rng);
	double cumulative = 0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		cumulative += p->birth;
		if( cumulative  > threshhold ){ 
			if(gsl_rng_uniform(rng) > MU){
				++p->popsize;
			} else {
				double step = gsl_ran_gaussian_ziggurat(rng,SIGMA_MU);
				pop pop_i = {1, p->trait+step, 0,0, p->trait, 0};
				poplist.push_back(pop_i);
			}
			break;
		}
		cumulative += p->death;
		if( cumulative  > threshhold ){ 
			--p->popsize;
			if( p->popsize == 0 ){
				poplist.erase( p ); 
			}
			break;
		}
	}
}


