qsub -t 1-211 -o /u/project/gandalm/cindywen/isoform_twas/star/log/job_out_walker_new_pass2 step-03-star-2ndpass.sh /u/project/gandalm/cindywen/isoform_twas/star/sample_list_walker_final211_vcf.txt walker

qsub -t 1-120 -o /u/project/gandalm/cindywen/isoform_twas/star/log/job_out_obrien_new_pass2 step-03-star-2ndpass.sh /u/project/gandalm/cindywen/isoform_twas/star/sample_list_obrien_vcf.txt obrien

qsub -t 1-116 -o /u/project/gandalm/cindywen/isoform_twas/star/log/job_out_werling_new_pass2 step-03-star-2ndpass.sh /u/project/gandalm/cindywen/isoform_twas/star/sample_list_werling_vcf.txt werling

qsub -t 1-44 -o /u/project/gandalm/cindywen/isoform_twas/star/log/job_out_libd_new_pass2 step-03-star-2ndpass.sh /u/project/gandalm/cindywen/isoform_twas/star/sample_list_libd_final44_vcf.txt libd

qsub -t 1-163 -o /u/project/gandalm/cindywen/isoform_twas/star/log/job_out_hdbr_new_pass2 step-03-star-2ndpass.sh /u/project/gandalm/cindywen/isoform_twas/star/sample_list_hdbr_final163_vcf.txt hdbr
