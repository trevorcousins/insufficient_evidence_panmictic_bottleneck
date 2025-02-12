for zpop in YRI ESN MSL GBR TSI CHB JPT
    do for trend_1 in 0 ;
        do for trend_2 in 1  
            do for ridge in 750 1000 5000
                do for folded in True False;
                    do
                    outprefix=/home/tc557/rds/hpc-work/mushi_1KGP_inference_241011/${zpop}_${trend_1}_${trend_2}_${ridge}_${folded}
                    outfileflag=${outprefix}SFS_expected_observed.png
                    if [ ! -f ${outfileflag} ]; then
                        echo file DOES NOT exist
                        commandline="python /home/tc557/falsifying_bottleneck/infer_mushi_model_241011.py -pop ${zpop} -out_prefix ${outprefix} -trend_1 ${trend_1} -trend_2 ${trend_2} -ridge ${ridge} -mu 1.25e-08 -pts 200 -most_ancient_gens 5e+04 -maxiter 300 -folded $folded"
                        echo $commandline
                        sbatch --wrap="$commandline" -p icelake -t 00:02:00 -A DURBIN-SL2-CPU --mem=2G
                    else
                        echo file DOES exist
                    fi
                    done
                done
            done
        done
    done
