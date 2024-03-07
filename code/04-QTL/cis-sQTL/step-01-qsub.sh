qsub -t 1-211 -o /u/project/gandalm/cindywen/isoform_twas/star/log/job_out_walker_new step-01-star-1stpass.sh /u/project/gandalm/cindywen/isoform_twas/star/sample_list_walker_final211.txt walker 

qsub -t 1-120 -o /u/project/gandalm/cindywen/isoform_twas/star/log/job_out_obrien_new step-01-star-1stpass.sh /u/project/gandalm/cindywen/isoform_twas/star/sample_list_obrien.txt obrien 

qsub -t 1-116 -o /u/project/gandalm/cindywen/isoform_twas/star/log/job_out_werling_new step-01-star-1stpass.sh /u/project/gandalm/cindywen/isoform_twas/star/sample_list_werling.txt werling 

qsub -t 1-44 -o /u/project/gandalm/cindywen/isoform_twas/star/log/job_out_libd_new step-01-star-1stpass.sh /u/project/gandalm/cindywen/isoform_twas/star/sample_list_libd_final44.txt libd 

qsub -t 1-163 -o /u/project/gandalm/cindywen/isoform_twas/star/log/job_out_hdbr_new step-01-star-1stpass.sh /u/project/gandalm/cindywen/isoform_twas/star/sample_list_hdbr_final163.txt hdbr 
