/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

#include "data_utility.h"
#include "irtree.h"

#include <unordered_map>

#ifndef WIN32

/*
 * GetCurTime is used to get the current running time in the current process.
 *
 * @Param curTime gets back the time information.
 *
 * @Return void.
 */
void GetCurTime(rusage *curTime) {
	int ret = getrusage(RUSAGE_SELF, curTime);
	if (ret != 0) {
		fprintf(stderr,
				"The running time info couldn't be collected successfully.\n");
		//FreeData( 2);
		exit(0);
	}
}

/*
 * GetTime is used to get the 'float' format time from the start and end rusage structure.
 *
 * @Param timeStart, timeEnd indicate the two time points.
 * @Param userTime, sysTime get back the time information.
 *
 * @Return void.
 */
void GetTime(struct rusage *timeStart, struct rusage *timeEnd, float *userTime,
		float *sysTime) {
	(*userTime) =
			((float) (timeEnd->ru_utime.tv_sec - timeStart->ru_utime.tv_sec))
					+ ((float) (timeEnd->ru_utime.tv_usec
							- timeStart->ru_utime.tv_usec)) * 1e-6;
	(*sysTime) =
			((float) (timeEnd->ru_stime.tv_sec - timeStart->ru_stime.tv_sec))
					+ ((float) (timeEnd->ru_stime.tv_usec
							- timeStart->ru_stime.tv_usec)) * 1e-6;
}

#endif

/*
 *	IRTree_read_config reads the configuration info fro constructing the IRTree.
 */
IRTree_config_t* read_config_irtree() {
	//char des[MAX_DESCRIPTION_LENG];
	FILE *c_fp;

	IRTree_config_t *cfg = (IRTree_config_t*) malloc(sizeof(IRTree_config_t));

	if ((c_fp = fopen(CONFIG_FILE, "r")) == NULL) {
		fprintf(stderr, "The config file cannot be opened.\n");
		exit(0);
	}

	//reads the configuration info.
	fscanf(c_fp, "%s%s%s", cfg->loc_file, cfg->doc_file, cfg->tree_file);
	fscanf(c_fp, "%i%i%i", &cfg->obj_n, &cfg->key_n, &cfg->dim); //data related.
	//fscanf( c_fp, "%i%i%i", &(cfg->M), &(cfg->m), &( cfg->split_opt));//R-tree related.
	fscanf(c_fp, "%i", &(cfg->split_opt)); //R-tree related.
	//fscanf( c_fp, "%i", &cfg->key_n);									//IF related.

	fclose(c_fp);

	return cfg;
}

/*
 *	Read the configuration for the Co-location mining problem.
 *  return an colocation_config_t pointer pointing to an object that storing the config
 */
colocation_config_t* read_config_colocation() {
	colocation_config_t *cfg;
	FILE *c_fp;

	cfg = (colocation_config_t*) malloc(sizeof(colocation_config_t));
	memset(cfg, 0, sizeof(colocation_config_t));

	if ((c_fp = fopen(COLOCATION_CONFIG_FILE, "r")) == NULL) {
		fprintf(stderr, "The colocation_config file cannot be opened.\n");
		exit(0);
	}

	//algorithm option.
	fscanf(c_fp, "%i%i", &cfg->alg_opt, &cfg->cost);

	//data.
	fscanf(c_fp, "%i%i%s", &cfg->dim, &cfg->obj_n, cfg->loc_file);

	fscanf(c_fp, "%i%s", &cfg->key_n, cfg->doc_file);

	fscanf(c_fp, "%i%s", &cfg->tree_tag, cfg->tree_file);

	fscanf(c_fp, "%i", &cfg->query_n);

	//colocation pattern mining.
	fscanf(c_fp, "%f", &cfg->min_sup);
	//fscanf(c_fp, "%f", &cfg->min_conf);
	fscanf(c_fp, "%f", &cfg->dist_thr);


	fclose(c_fp);

	return cfg;
}

/*
 *	Add a key to the keyword list. //!!k_node_v is the pointer pointing the last element of the list!!
 */
void add_keyword_entry(k_node_t *&k_node_v, KEY_TYPE key) {
	k_node_v->next = (k_node_t*) malloc(sizeof(k_node_t));
	memset(k_node_v->next, 0, sizeof(k_node_t));
	k_node_v->next->key = key;
	k_node_v = k_node_v->next;

	/*s*/
	stat_v.memory_v += sizeof(k_node_t);
	if (stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
}

/*
 *	Copy the k_list info of @k_head2 to @k_head1.
 */
void copy_k_list(k_node_t *k_head1, k_node_t *k_head2) {
	k_node_t *k_node_iter1, *k_node_iter2;

	k_node_iter1 = k_head1;
	k_node_iter2 = k_head2->next;
	while (k_node_iter2 != NULL) {
		add_keyword_entry(k_node_iter1, k_node_iter2->key);

		k_node_iter2 = k_node_iter2->next;
	}
}

/*
 *	Print the a list of keywords @k_head in @o_fp.
 */
void print_k_list(k_node_t *k_head, FILE *o_fp) {
	k_node_t *k_node_iter;

	k_node_iter = k_head->next;
	while (k_node_iter != NULL) {
		fprintf(o_fp, "%.0lf  ", k_node_iter->key);

		k_node_iter = k_node_iter->next;
	}

	fprintf(o_fp, "\n");
}

/****
 *	Print the location @k_head in @o_fp.
 ****/
void print_loc(loc_t *loc_v, FILE *o_fp) {
	int i;
	for (i = 0; i < loc_v->dim; i++) {
		fprintf(o_fp, "%0.4lf  ", loc_v->coord[i]);
	}
	fprintf(o_fp, "\n");
}

/*
 *	Print the statistics maintained in stat_v.
 */

void print_colocation_stat(colocation_config_t *cfg, int cnt) {
	FILE *s_fp;

	if (!(s_fp = fopen(COLOCATION_STAT_FILE, "w"))) {
		fprintf(stderr, "Cannot open the coskq_stat file.\n");
		exit(0);
	}

	//	//average cost.
	//	fprintf( s_fp, "%lf\n", stat_v.aver_cost);
	//	//average size.
	//	fprintf( s_fp, "%lf\n\n", stat_v.aver_size);

	//time.
	fprintf(s_fp, "%f\n%f\n\n", stat_v.irtree_build_time, stat_v.q_time+ stat_v.pre_time);

	//memory.
	fprintf(s_fp, "%f\n", stat_v.memory_max / (1024 * 1024));

	//IR-tree memory.
	fprintf(s_fp, "%f\n\n", stat_v.tree_memory_max / (1024 * 1024));

	//Method 4 related
	fprintf(s_fp, "%lf\n%lf\n%lf\n%lf\n%lf\n\n", stat_v.S0_sum,
			stat_v.S3_sum,  stat_v.S1_sum, stat_v.S2_sum,
			stat_v.S5_sum);

	fprintf(s_fp, "%lf\n%lf\n%lf\n%lf\n%lf\n", stat_v.S0_time,
			stat_v.S3_time, stat_v.S1_time, stat_v.S2_time,
			stat_v.S5_time);

	fclose(s_fp);
}

/*
 *	Allocate the memory for an object.
 */
void alloc_obj(obj_t *obj_v, int dim) {
	//obj_v = ( obj_t*)malloc( sizeof( obj_t));
	//memset( obj_v, 0, sizeof( obj_t));

	obj_v->MBR = (range*) malloc(dim * sizeof(range));
	memset(obj_v->MBR, 0, dim * sizeof(range));

	//    obj_v->N_o_f = new std::unordered_map<FEA_TYPE, int>();
	//
	//    obj_v->N_o_f2 = new std::unordered_map<FEA_TYPE, int>();
	//
	//    obj_v->frac_f = new std::unordered_map<FEA_TYPE, float>();
	//
}

/*
 *	Read the data based on the IRTree_config_t info.
 */
data_t* read_data_irtree(IRTree_config_t *cfg) {
	int i, j;
	KEY_TYPE key;
	char des;
	char keys[TEXT_COL_MAX];
	char *tok;

	k_node_t *k_node_v;
	FILE *i_fp;

	data_t *data_v = (data_t*) malloc(sizeof(data_t));
	memset(data_v, 0, sizeof(data_t));

	data_v->dim = cfg->dim;
	data_v->obj_n = cfg->obj_n;
	data_v->key_n = cfg->key_n;

	data_v->obj_v = (obj_t*) malloc(sizeof(obj_t) * data_v->obj_n);
	memset(data_v->obj_v, 0, sizeof(obj_t) * data_v->obj_n);

	//data_v->key_freq_v = bst_ini( );

	//Read the loc info.
	if ((i_fp = fopen(cfg->loc_file, "r")) == NULL) {
		fprintf(stderr, "Cannot open the loc file.\n");
		exit(0);
	}

	for (i = 0; i < data_v->obj_n; i++) {
		alloc_obj(data_v->obj_v + i, data_v->dim);

		fscanf(i_fp, "%i", &(data_v->obj_v[i].id));

		for (j = 0; j < data_v->dim; j++) {
			fscanf(i_fp, "%c%f", &des, &(data_v->obj_v[i].MBR[j].min));
			data_v->obj_v[i].MBR[j].max = data_v->obj_v[i].MBR[j].min;
		}
	}

	fclose(i_fp);

	//Read the keywords info.
	if ((i_fp = fopen(cfg->doc_file, "r")) == NULL) {
		fprintf(stderr, "Cannot open the doc file.\n");
		exit(0);
	}

	for (i = 0; i < data_v->obj_n; i++) {
		//        k_node_v = ( k_node_t*)malloc( sizeof( k_node_t));
		//        memset( k_node_v, 0, sizeof( k_node_t));

		//        data_v->obj_v[ i].k_head = k_node_v;

		fgets(keys, TEXT_COL_MAX, i_fp);

		tok = strtok(keys, " ,");
		while ((tok = strtok(NULL, " ,"))) {
			key = atoi(tok);

			data_v->obj_v[i].fea = key;
			//    add_keyword_entry( k_node_v, key);
			/*
			 //Update the frequency info of the keyword.
			 bst_node_v = bst_search( data_v->key_freq_v, key);
			 if( bst_node_v == NULL)
			 {
			 bst_node_v = ( bst_node_t*)malloc( sizeof( bst_node_t));
			 memset( bst_node_v, 0, sizeof( bst_node_t));
			 
			 bst_node_v->key = key;
			 bst_node_v->freq = 1.0 / data_v->obj_n;
			 
			 bst_insert( data_v->key_freq_v, bst_node_v);
			 }
			 else
			 bst_node_v->freq += 1.0 / data_v->obj_n;
			 */
		}
	}

	fclose(i_fp);

	return data_v;
}

/*
 *
 */
data_t* alloc_data(int num) {
	data_t *data_v;

	data_v = (data_t*) malloc(sizeof(data_t));
	memset(data_v, 0, sizeof(data_t));

	data_v->obj_n = num;

	data_v->obj_v = (obj_t*) malloc(sizeof(obj_t) * data_v->obj_n);
	memset(data_v->obj_v, 0, sizeof(obj_t) * data_v->obj_n);

	return data_v;
}


/*
 *	Read the data based on the colocation_config_t info.
 */
data_t* read_data_colocation(colocation_config_t *cfg) {
	int i, j;
	KEY_TYPE key;
	char des;
	char keys[TEXT_COL_MAX];
	char *tok;
	data_t *data_v;

	FILE *i_fp;

	data_v = alloc_data(cfg->obj_n);
	// performed inside alloc_data function: data_v->obj_n = cfg->obj_n;

	data_v->dim = cfg->dim;
	data_v->key_n = cfg->key_n;

	//data_v->key_freq_v = bst_ini( );

	///--------------------------------------------------
	//Read the loc info.
	/// from the loc file
	if ((i_fp = fopen(cfg->loc_file, "r")) == NULL) {
		fprintf(stderr, "Cannot open the loc file.\n");
		exit(0);
	}

//    for (i = 0; i < data_v->obj_n; i++) {
//        alloc_obj(data_v->obj_v + i, data_v->dim);
//
//        fscanf(i_fp, "%i", &(data_v->obj_v[i].id));
//        if(data_v->obj_v[i].id<10000)
//        printf("%d\n",data_v->obj_v[i].id);
//
//        for (j = 0; j < data_v->dim; j++) {
//            fscanf(i_fp, "%c%f", &des, &(data_v->obj_v[i].MBR[j].min));
//            data_v->obj_v[i].MBR[j].max = data_v->obj_v[i].MBR[j].min;
//        }
//    }
	for (i = 0; i < data_v->obj_n; i++) {

		alloc_obj(data_v->obj_v + i, data_v->dim);

		fgets(keys, TEXT_COL_MAX, i_fp);
//            printf("keys: %s\n",keys);
		tok = strtok(keys, ",");
		data_v->obj_v[i].id = atoi(tok);
		int j = 0;
		while ((tok = strtok(NULL, ","))) {

			key = atof(tok);
			data_v->obj_v[i].MBR[j].min = key;
			data_v->obj_v[i].MBR[j].max = data_v->obj_v[i].MBR[j].min;
			j++;
		}
	}
	fclose(i_fp);
	///--------------------------------------------------
	//Read the keywords info.
	///from the doc file
	if ((i_fp = fopen(cfg->doc_file, "r")) == NULL) {
		fprintf(stderr, "Cannot open the doc file.\n");
		exit(0);
	}

	for (i = 0; i < data_v->obj_n; i++) {
		fgets(keys, TEXT_COL_MAX, i_fp);

		tok = strtok(keys, ",");
//        while ((tok = strtok(NULL, ","))) {
//            key = atoi(tok);
//            data_v->obj_v[i].fea = key;
//        }
		//fea
		tok = strtok(NULL, ",");
		data_v->obj_v[i].fea = atoi(tok);
		//weight
		tok = strtok(NULL, ",");
		if (tok != NULL)
			data_v->obj_v[i].weight = atof(tok);
		else
			data_v->obj_v[i].weight = 1;
//			printf("data_v->obj_v[i].weight: %.2lf\n",data_v->obj_v[i].weight );

	}
	fclose(i_fp);

	return data_v;
}

/*
 *	Release a list.
 */
void release_k_list(k_node_t *k_node_v) {
	k_node_t *k_node_v1;
	while (k_node_v->next != NULL) {
		k_node_v1 = k_node_v;
		k_node_v = k_node_v->next;
		free(k_node_v1);
	}
	free(k_node_v);
}

/*
 *	IRTree_free_data release the data read by 'IRTree_read_data'.
 */
void release_data(data_t *data_v) {
	int i;

	for (i = 0; i < data_v->obj_n; i++) {
		free(data_v->obj_v[i].MBR);
	}

	free(data_v->obj_v);
	free(data_v);
}

/*
 *	Generates a number between [min, max] randomly.
 */
int rand_i(int min, int max) {
	int i_rand;
	float ratio;

	if (min > max)
		return 0;

	if (max - min > RAND_MAX) {
		//printf( "rand_i: [%i, %i] out of range.\n", min, max);
		//exit( 0);

		ratio = rand_f(0, 1);
		i_rand = min + (int) (ratio * (max - min));

		return i_rand;
	}

	//srand( time( NULL));

	if (min == max)
		i_rand = 0;
	else
		i_rand = rand() % (max - min + 1);

	return i_rand + min;
}

/*
 *	Generate a random float value in [min_v, max_v].
 */
float rand_f(float min_v, float max_v) {
	float rand_v;

	if (min_v > max_v)
		return 0;

	//srand( time(NULL));

	rand_v = float(rand()) / RAND_MAX; //rand_v falls in [0,1] randomly.

	return (max_v - min_v) * rand_v + min_v;
}

/*
 *	Generate a float value that satisfy the normal/gaussian distribution N[mean, s_var^2].
 *
 *	Polar form method, modified from the Internet.
 */
float gaussian_f(float mean, float s_var) {
	float x1, x2, w, y1, y2;

	do {
		x1 = 2 * rand_f(0, 1) - 1;
		x2 = 2 * rand_f(0, 1) - 1;
		w = x1 * x1 + x2 * x2;
	} while (w >= 1);

	w = float(sqrt(float(-2.0 * log(w)) / w));
	y1 = x1 * w; //y1 follows N[0, 1^2].

	y2 = y1 * s_var + mean; //y2 follows N[mean, s_var^2].

	return y2;
}

/*
 *	Check whether an variable has been generated before.
 *
 */
bool is_old(int *rand_v, int cur, int var) {
	int i;

	for (i = 0; i < cur; i++) {
		if (rand_v[i] == var)
			return true;
	}

	return false;
}

/*
 *
 */
range* collect_data_range(data_t *data_v) {
	int i, j;
	range *MBR;

	MBR = (range*) malloc(data_v->dim * sizeof(range));
	memset(MBR, 0, data_v->dim * sizeof(range));

	for (i = 0; i < data_v->dim; i++) {
		MBR[i].min = FLT_MAX;
		MBR[i].max = -FLT_MAX;

		for (j = 0; j < data_v->obj_n; j++) {
			if (data_v->obj_v[j].MBR[i].min < MBR[i].min)
				MBR[i].min = data_v->obj_v[j].MBR[i].min;

			if (data_v->obj_v[j].MBR[i].max > MBR[i].max)
				MBR[i].max = data_v->obj_v[j].MBR[i].max;
		}
	}

	return MBR;
}

/*
 *
 */
void print_data(data_t *data_v, syn_config_t *s_cfg) {
	int i, j;

	FILE *loc_fp, *doc_fp;

	if ((loc_fp = fopen(s_cfg->loc_file, "w")) == NULL
			|| (doc_fp = fopen(s_cfg->doc_file, "w")) == NULL) {
		printf("Cannot open loc/doc files.\n");
		exit(0);
	}

	//Print the loc info.
	for (i = 0; i < data_v->obj_n; i++) {
		fprintf(loc_fp, "%i", data_v->obj_v[i].id);
		for (j = 0; j < data_v->dim; j++)
			fprintf(loc_fp, ",%f", data_v->obj_v[i].MBR[j].min);

		fprintf(loc_fp, "\n");
	}

	//Print the doc info.
	for (i = 0; i < data_v->obj_n; i++) {
		fprintf(doc_fp, "%i", data_v->obj_v[i].id);
		fprintf(doc_fp, "%i", (int) data_v->obj_v[i].fea);
		fprintf(doc_fp, "\n");
	}

	fclose(loc_fp);
	fclose(doc_fp);
}

//Following the paper "Discovering Spatial Co-location Patterns" in SSTD01
void gen_syn_data(syn_config_t *s_cfg) {
	//parameters
	int Ncoloc = s_cfg->Ncoloc; //number of colocation

	int lambda1 = s_cfg->lambda1; //size of each colocation pattern
	int lambda2 = s_cfg->lambda2; //size of each table instance of a colocation pattern

	int seed1 = s_cfg->seed1; //random seed for lambda1
	int seed2 = s_cfg->seed2; //random seed for lambda2

	double D1 = s_cfg->D1, D2 = s_cfg->D2;
	double d = s_cfg->d;

	double r_noise_fea = s_cfg->r_noise_fea; //number of *additional* features for noise / number of original features
	double r_global = s_cfg->r_global; //number of *additional* instances from noise features
	double r_local = s_cfg->r_local; //number of *additional* instances from original features

	// printf("D1:%f\tD2:%f\n",D1,D2);
	//----------------

	int objectID = 1; //starting objectID
	int feaID = 1; //starting feature ID

	int num_instance = 0; //number of instances
	FILE *loc_fp, *doc_fp;

	if ((loc_fp = fopen(s_cfg->loc_file, "w")) == NULL
			|| (doc_fp = fopen(s_cfg->doc_file, "w")) == NULL) {
		printf("Cannot open loc/doc files.\n");
		exit(0);
	}

	std::default_random_engine generator1, generator2;
	std::poisson_distribution<int> distribution1(lambda1), distribution2(
			lambda2);
	generator1.seed(seed1); //[1,infty)
	generator2.seed(seed2); //[1,infty)

	int *patternSize = new int[Ncoloc];

	//generate N patterns size
	for (int i = 0; i < Ncoloc; i++) {
		patternSize[i] = distribution1(generator1);

		int tableSize = distribution2(generator2);
		//        printf("tableSize:%d\tpatternSize:%d\n",tableSize, patternSize[i]);

		num_instance += tableSize * patternSize[i];

		//create objects for each row instance in the table
		for (int j = 0; j < tableSize; j++) {
			double rectX = rand_f(0, D1 - d); //bottom left
			double rectY = rand_f(0, D2 - d);

			//each object correspond to one feature
			for (int k = 0; k < patternSize[i]; k++) {
				int fea = feaID + k;
				double locX = rectX + rand_f(0, d);
				double locY = rectY + rand_f(0, d);
				//                printf("%d\t%d\t%0.3lf\t%0.3lf\n",objectID,fea,locX, locY);
				fprintf(loc_fp, "%d,%f,%f\n", objectID, locX, locY);
				fprintf(doc_fp, "%d,%d\n", objectID, fea);
				objectID++;
			}
			//            printf("\n");
		}
		feaID += patternSize[i];

		//        printf("===============\n");
	}

	//    printf("objID:%d\tnum_instance:%d\tfeaID:%d\n", objectID, num_instance, feaID);

	//    printf("------------ local noise -----------\n");

	//generate local noise
	for (int i = 0; i < r_local * num_instance; i++) {
		//pick from existing features
		int fea = rand_i(1, feaID);
		double locX = rand_f(0, D1);
		double locY = rand_f(0, D2);
		//        printf("%d\t%d\t%0.3lf\t%0.3lf\n",objectID,fea,locX, locY);
		fprintf(loc_fp, "%d,%f,%f\n", objectID, locX, locY);
		fprintf(doc_fp, "%d,%d\n", objectID, fea);
		objectID++;
	}

	//    printf("------------ global noise -----------\n");

	int originalFeaID = feaID;
	feaID = int(double(feaID) * (1.0 + r_noise_fea));

	//generate global noise
	for (int i = 0; i < r_global * num_instance; i++) {
		//pick from new features
		int fea = rand_i(originalFeaID, feaID);
		double locX = rand_f(0, D1);
		double locY = rand_f(0, D2);
		//   printf("%d\t%d\t%0.3lf\t%0.3lf\n",objectID,fea,locX, locY);
		fprintf(loc_fp, "%d,%f,%f\n", objectID, locX, locY);
		fprintf(doc_fp, "%d,%d\n", objectID, fea);
		objectID++;
	}

	printf("total number of obj:%d\ntotal number of feature:%d\n", --objectID,
			feaID);
	printf("*fea[%d,%d] are original, [%d,%d] are noise\n", 1, originalFeaID,
			originalFeaID + 1, feaID);

	fclose(loc_fp);
	fclose(doc_fp);
}

/*
 *
 */
void batch_gen_syn_data() {
	int ins_n;
	syn_config_t s_cfg;
	FILE *c_fp;

	if ((c_fp = fopen(SYN_CONFIG_FILE, "r")) == NULL) {
		printf("Cannot open the syn_config file.\n");
		exit(0);
	}

	ins_n = 1;
	while (fscanf(c_fp, "%s%s", s_cfg.loc_file, s_cfg.doc_file) != EOF) {
		printf("#Instance %i ...\n", ins_n++);

		//read the config.
		fscanf(c_fp, "%i%i%i", &s_cfg.Ncoloc, &s_cfg.lambda1, &s_cfg.lambda2);
		fscanf(c_fp, "%i%i", &s_cfg.seed1, &s_cfg.seed2);
		fscanf(c_fp, "%lf%lf%lf", &s_cfg.D1, &s_cfg.D2, &s_cfg.d);

		fscanf(c_fp, "%lf%lf%lf", &s_cfg.r_noise_fea, &s_cfg.r_global,
				&s_cfg.r_local);

		fscanf(c_fp, "%lf%lf", &s_cfg.m_overlap, &s_cfg.m_clump);

		//generate the synthetic data.
		gen_syn_data(&s_cfg);
	}
}

//----------------------------------

/*
 * generate synethetic datasets in batch, each row in syn_config.txt correposnd to a dataset
 */
void batch_gen_syn_data2() {
	int ins_n;
	syn_config_t s_cfg;
	FILE *c_fp;

	if ((c_fp = fopen(SYN_CONFIG_FILE, "r")) == NULL) {
		printf("Cannot open the syn_config file.\n");
		exit(0);
	}

	ins_n = 1;
	while (fscanf(c_fp, "%s%s", s_cfg.loc_file, s_cfg.doc_file) != EOF) {
		printf("#Instance %i ...\n", ins_n++);

		//read the config.
		fscanf(c_fp, "%i%i%i", &s_cfg.Ncoloc, &s_cfg.lambda1, &s_cfg.lambda2);
		fscanf(c_fp, "%i%i", &s_cfg.seed1, &s_cfg.seed2);
		fscanf(c_fp, "%lf%lf%lf", &s_cfg.D1, &s_cfg.D2, &s_cfg.d);

		fscanf(c_fp, "%lf%lf%lf", &s_cfg.r_noise_fea, &s_cfg.r_global,
				&s_cfg.r_local);

		fscanf(c_fp, "%lf%lf%i", &s_cfg.m_overlap, &s_cfg.m_clump,
				&s_cfg.m_clump_type);

		//generate the synthetic data.
		gen_syn_data2(&s_cfg);
	}
}

//Following the paper "Discovering colocation patterns from spatial data sets: A general approach" in TKDE04
void gen_syn_data2(syn_config_t *s_cfg) {
	//parameters
	int Ncoloc = s_cfg->Ncoloc; //number of colocation
	int lambda1 = s_cfg->lambda1; //size of each colocation pattern
	int lambda2 = s_cfg->lambda2; //size of each table instance of a colocation pattern
	int seed1 = s_cfg->seed1; //random seed for lambda1
	int seed2 = s_cfg->seed2; //random seed for lambda2
	double D1 = s_cfg->D1; //D2 should always same as D1 and thus D2 is not used
	double d = s_cfg->d;
	double r_noise_fea = s_cfg->r_noise_fea; //number of *additional* features for noise / number of original features
	double r_local = s_cfg->r_local;
	double r_global = s_cfg->r_global; //number of *additional* instances from noise features
	double m_overlap = s_cfg->m_overlap;
	double m_clump = s_cfg->m_clump;
	int m_clump_type = s_cfg->m_clump_type;

	//----------------
	int objectID = 1; //starting objectID
	int feaID = 1; //starting feature ID
	int num_instance = 0; //number of instances

	FILE *loc_fp, *doc_fp;
	if ((loc_fp = fopen(s_cfg->loc_file, "w")) == NULL
			|| (doc_fp = fopen(s_cfg->doc_file, "w")) == NULL) {
		printf("Cannot open loc/doc files.\n");
		exit(0);
	}

	std::default_random_engine generator1, generator2, generator3;
	std::poisson_distribution<int> distribution1(lambda1), distribution2(
			lambda2), distribution3(m_clump);
	generator1.seed(seed1); //[1,infty)
	generator2.seed(seed2); //[1,infty)
	generator3.seed(seed2);

	int *patternSize = new int[Ncoloc]; //storing the size |C| (number of features) of each pattern

	//generate N patterns size
	for (int i = 0; i < Ncoloc; i++) {
		patternSize[i] = distribution1(generator1);

		//generate m_overlap maximal colocation for each core colocation
		for (int a = 0; a < m_overlap; a++) {
			int tableSize = distribution2(generator2);
			printf("%d\t%d\t\t%d\n", patternSize[i] + 1, tableSize, feaID);

			int extrafeaID = feaID + patternSize[i] + a; //the one more spatial feature in each maximal colocation

			//---
			//generate m_clump value for each feature in this pattern
			//m_clump_type ==0: all same fxied value m_clump
			//1: random
			int *m_clump_value = new int[patternSize[i] + 1];
			for (int k = 0; k < patternSize[i] + 1; k++) {
				if (m_clump_type == 0)
					m_clump_value[k] = m_clump;
				else if (m_clump_type == 1)
					m_clump_value[k] = rand_i(1, m_clump);
				printf("%d\t", m_clump_value[k]);
			}
			printf("\n");
			//---

			//create objects for each row instance in the table
			for (int j = 0; j < tableSize; j++) {
				double rectX = rand_f(0, D1 - d); //bottom left
				double rectY = rand_f(0, D1 - d);

				//each object is associated with one feature
				for (int k = 0; k < patternSize[i]; k++) {
					for (int b = 0; b < m_clump_value[k]; b++) {
						gen_syn_obj(loc_fp, doc_fp, objectID, feaID + k, rectX,
								rectY, d);
						//generate additional objects right next to the current cell
						//						if(m_clump_type==1)
						//							gen_syn_obj(loc_fp, doc_fp, objectID, feaID + k, rectX+d, rectY, d);
					}
				}
				//create an object for the one more spatial feature
				for (int b = 0; b < m_clump_value[patternSize[i]]; b++) {
					gen_syn_obj(loc_fp, doc_fp, objectID, extrafeaID, rectX,
							rectY, d);
					//					if(m_clump_type==1)
					//						gen_syn_obj(loc_fp, doc_fp, objectID, extrafeaID, rectX+d, rectY, d);
				}
			}
			num_instance += tableSize * (patternSize[i] + 1) * m_clump;
		}
		feaID += patternSize[i] + m_overlap;
		//        printf("\n");
	}
	//	if(m_clump_type==1)
	//		num_instance *=2;

	//    printf("objID:%d\tnum_instance:%d\tfeaID:%d\n", objectID, num_instance, feaID);
	//------generate local noise-------
	gen_syn_noise(r_local * num_instance, loc_fp, doc_fp, objectID, 1, feaID, 0,
			0, D1);

	//    printf("------------ global noise -----------\n");
	int originalFeaID = feaID;
	feaID = int(double(feaID) * (1.0 + r_noise_fea));
	gen_syn_noise(r_global * num_instance, loc_fp, doc_fp, objectID,
			originalFeaID, feaID, 0, 0, D1);

	printf("total number of obj:%d\ntotal number of feature:%d\n", --objectID,
			feaID);
	printf("*fea[%d,%d] are original, [%d,%d] are noise\n", 1,
			originalFeaID - 1, originalFeaID, feaID);

	fclose(loc_fp);
	fclose(doc_fp);
}

/**
 * generate *one* object id @objectID with feature @fea randomly in [rectX, rectY] to [rectX+d, rectY+d]
 */
void gen_syn_obj(FILE *loc_fp, FILE *doc_fp, int &objectID, int fea,
		double rectX, double rectY, double d) {
	double locX = rectX + rand_f(0, d);
	double locY = rectY + rand_f(0, d);
	//                printf("%d\t%d\t%0.3lf\t%0.3lf\n",objectID,fea,locX, locY);
	fprintf(loc_fp, "%d,%f,%f\n", objectID, locX, locY);
	fprintf(doc_fp, "%d,%d\n", objectID, fea);
	objectID++;
}

/**
 * generate local noise and global noise
 * @num_noise_obj noise objects
 * each randomly pick fea in range [feaStart, feaEnd]
 *
 */
void gen_syn_noise(int num_noise_obj, FILE *loc_fp, FILE *doc_fp, int &objectID,
		int feaStart, int feaEnd, double rectX, double rectY, double d) {
	if (num_noise_obj > 0) {
		for (int i = 0; i < num_noise_obj; i++) {
			//pick from new features
			int fea = rand_i(feaStart, feaEnd);
			gen_syn_obj(loc_fp, doc_fp, objectID, fea, 0, 0, d);
		}
	}
}

/*
 * generate scalability test dataset from real
 *
 */

void gen_scalability_data(data_t *data_v) {
	FILE *loc_fp, *doc_fp;
	if ((loc_fp = fopen("sca_loc.txt", "w")) == NULL
			|| (doc_fp = fopen("sca_doc.txt", "w")) == NULL) {
		printf("Cannot open loc/doc files.\n");
		exit(0);
	}

	//    int total = 1000000;
	//    for (int cnt = 1; cnt <= total; cnt++) {
	//pick an object and copy its feature
	//        int i = rand_i(0, data_v->obj_n - 1);

	int cnt = 1;
	for (int i = 0; i < data_v->obj_n; i++) {
		//number of times
		for (int t = 0; t < 28; t++) {

			//pick two objects and use mid point as new location
			int j = rand_i(0, data_v->obj_n - 1);
//                int k = rand_i(0, data_v->obj_n - 1);

			int fea = data_v->obj_v[i].fea;
			double x = (data_v->obj_v[j].MBR[0].min
					+ data_v->obj_v[i].MBR[0].min) / 2.0;
			double y = (data_v->obj_v[j].MBR[1].min
					+ data_v->obj_v[i].MBR[1].min) / 2.0;

			//print a new object
			fprintf(loc_fp, "%d,%f,%f\n", cnt, x, y);
			fprintf(doc_fp, "%d,%d\n", cnt, fea);
			cnt++;
		}
	}
	fclose(loc_fp);
	fclose(doc_fp);
}

/* scalability based on number of features
 * for each object, we randomize its feature
 */
void gen_scalability_data2(data_t *data_v) {
	FILE *loc_fp, *doc_fp;
	if ((loc_fp = fopen("sca_loc.txt", "w")) == NULL
			|| (doc_fp = fopen("sca_doc.txt", "w")) == NULL) {
		printf("Cannot open loc/doc files.\n");
		exit(0);
	}

	int cnt = 1;
	for (int i = 0; i < data_v->obj_n; i++) {
		//pick an object and copy its feature
		int r = rand_i(0, 1) * 36;
		//pick two objects and use mid point as new location
//		int j = rand_i(0, data_v->obj_n - 1);
//		int k = rand_i(0, data_v->obj_n - 1);

		int fea = data_v->obj_v[i].fea + r;
//		double x = (data_v->obj_v[j].MBR[0].min +data_v->obj_v[k].MBR[0].min)/2.0;
//		double y = (data_v->obj_v[j].MBR[1].min +data_v->obj_v[k].MBR[1].min)/2.0;
		double x = data_v->obj_v[i].MBR[0].min;
		double y = data_v->obj_v[i].MBR[1].min;
		//print a new object
		fprintf(loc_fp, "%d,%f,%f\n", cnt, x, y);
		fprintf(doc_fp, "%d,%d\n", cnt, fea);
		cnt++;
	}
	fclose(loc_fp);
	fclose(doc_fp);

}
