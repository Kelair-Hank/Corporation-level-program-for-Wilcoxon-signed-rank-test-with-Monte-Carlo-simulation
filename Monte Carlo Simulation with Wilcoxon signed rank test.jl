# mean_correct_rates_of_alignment.csv: The file records mean correct_rate of each subread.

using HypothesisTests, CSV, DataFrames, Pipe, Statistics

significant_percent(v)=@pipe Statistics.mean(v .< 0.05)*100 |> round(_ , digits=2)

# original_dir: The directory where the original algorithm stores its mean correct_rate of each subread.
# test_dirs: The directories where the test algorithms store their mean correct_rate of each subread.
# mean_correct_rates_file: The file records mean correct_rate of each subread.
# test_algorithm_names: Tested algorithms names.
function monte_carlo_simulation_with_wilcox_signed_rank_test(original_dir::String,test_dirs::Vector{String},mean_correct_rates_file::String,test_algorithm_names::Vector{Symbol})
    wd1=pwd()
    cd(original_dir)
    mean_correct_rates=CSV.read(mean_correct_rates_file,DataFrame,copycols=true)
    correct_rates=mean_correct_rates.correct_rates

    for dir1=test_dirs
        cd(dir1)
        mean_correct_rates1=CSV.read(mean_correct_rates_file,DataFrame,copycols=true)
        correct_rates1=mean_correct_rates1.correct_rates
        correct_rates=[correct_rates correct_rates1]
    end

    correct_rates=DataFrames.DataFrame(correct_rates)
    DataFrames.dropmissing!(correct_rates)
    
    p_values=[Inf Inf Inf Inf]
    Wilcoxon_stats=[Inf Inf Inf Inf]
    for i=1:1000
        sample1=correct_rates[rand(1:nrow(correct_rates),50),:]

        sample_p_values=Float64[Inf,]
        sample_Wilcoxon_stats=Float64[Inf,]
        for j=2:ncol(correct_rates)
            test_result=SignedRankTest(sample1[:,1],sample1[:,j])
            sample_p_values=[sample_p_values pvalue(test_result)]
            sample_Wilcoxon_stats=[sample_Wilcoxon_stats test_result.W]
        end

        p_values=[p_values;sample_p_values]
        Wilcoxon_stats=[Wilcoxon_stats;sample_Wilcoxon_stats]
    end

    p_values=DataFrames.DataFrame(p_values)
    Wilcoxon_stats=DataFrames.DataFrame(Wilcoxon_stats)

    DataFrames.deleterows!(p_values,1)
    DataFrames.delete!(p_values,1)
    DataFrames.deleterows!(Wilcoxon_stats,1)
    DataFrames.delete!(Wilcoxon_stats,1)
    DataFrames.rename!(p_values,test_algorithm_names)
    DataFrames.rename!(Wilcoxon_stats,test_algorithm_names)

    p_values_mean=DataFrames.describe(p_values,:mean)
    Wilcoxon_stats_mean=DataFrames.describe(Wilcoxon_stats,:mean)
    
    p_values[!,:group]=repeat(["ab",],nrow(p_values))
    p_values_group=DataFrames.groupby(p_values,:group)
    p_values_significant=DataFrames.combine(p_values_group,test_algorithm_names .=> significant_percent)

    cd(wd1)
    DataFrames.delete!(p_values,ncol(p_values))
    DataFrames.delete!(p_values_significant,1)
    CSV.write("p_values.csv",p_values)
    CSV.write("Wilcoxon_stats.csv",Wilcoxon_stats)
    CSV.write("p_values_mean.csv",p_values_mean)
    CSV.write("Wilcoxon_stats_mean.csv",Wilcoxon_stats_mean)
    CSV.write("p_values_significant.csv",p_values_significant)
end

# The following directories may need to be changed each time. 
dir1="/mnt/projects/01.ccs/results/2/baseCall-27k.rar/2022-05-17/ccs/temp"
dir2="/mnt/projects/01.ccs/results/2/baseCall-27k.rar/2022-05-16/ccs/temp"
dir3="/mnt/projects/01.ccs/results/2/baseCall-27k.rar/2022-05-19/ccs/temp"
dir4="/mnt/projects/01.ccs/results/2/baseCall-27k.rar/2022-05-24/ccs/temp"

# The following file name and algorithm names may need to be changed each time.
monte_carlo_simulation_with_wilcox_signed_rank_test(dir1,String[dir2,dir3,dir4],"mean_correct_rates_of_alignment.csv",Symbol[:FAMSA,:Mafft,:MUSCLE])
