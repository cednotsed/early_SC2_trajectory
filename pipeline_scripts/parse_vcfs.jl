using CodecZlib


function parse_vcfs(in_path, out_path)
    # Open output file
    open(out_path, "w") do out_file
        # Headers
        println(out_file, "pos\tindel_flag\tref\talt\tqual\tmq\ttotal_depth\tad\tn_A\tn_T\tn_G\tn_C\tp_A\tp_T\tp_G\tp_C")

        # Open and read the VCF file
        open(in_path, "r") do in_file
            decoder = GzipDecompressorStream(in_file)

            for line in eachline(decoder)
                # Skip header lines
                if startswith(line, "#")
                    continue
                end
        
                # Split the VCF record line into fields
                fields = split(line, '\t')
                #@show fields 
                # Extract necessary information
                pos = fields[2]
                ref = fields[4]
                alts = fields[5]
                qual = fields[6]
                allele_names = vcat(ref, split(alts, ','))
                info_string = fields[8]
                indel_flag = split(info_string, ';')[1] == "INDEL"
                mq = split(info_string, "MQ=")[2]
                format_string = split(fields[10], ':')
                dp = parse(Int64, format_string[3])
                ad_string = format_string[4]
                ad = [parse(Int64, ss) for ss in split(ad_string, ',')]
        
                #@show pos
                #@show ref
                #@show alts
                #@show allele_names
                #@show ad
                
                # Initialize an array to aggregate AD counts for each of the 4 alleles (A, C, G, T)
                allele_counts = Dict('A' => 0, 'T' => 0, 'G' => 0, 'C' => 0)
                allele_freqs = Dict('A' => 0.0, 'T' => 0.0, 'G' => 0.0, 'C' => 0.0)

                # Assign allele counts
                if dp > 0
                    for i in 1:length(allele_names)
                        allele_counts[allele_names[i][1]] = ad[i]
                        allele_freqs[allele_names[i][1]] = ad[i] / dp
                    end
                end
        
                #@show allele_counts
                
                # Parse result line
                parsed_row = "$pos\t$indel_flag\t$ref\t$alts\t$qual\t$mq\t$dp\t$ad_string"
        
                # Add allele counts
                for key in alleles
                    allele_count = allele_counts[key]
                    parsed_row = parsed_row * "\t$allele_count"  
                end
                
                for key in alleles
                    allele_freq = allele_freqs[key]
                    parsed_row = parsed_row * "\t$allele_freq"
                end
        
                #@show parsed_row
        
                # Write lines
                println(out_file, parsed_row)
            end
        end
    end
end


# MAIN
in_dir = "/mnt/c/git_repos/wuhu_rooting/results/pipeline_out/vcf_out/"
out_dir = "/mnt/c/git_repos/wuhu_rooting/results/pipeline_out/parsed_vcfs/"

# Define possible alleles
global alleles = ['A', 'T', 'G', 'C']
 
# Multithread across multiple VCFs
file_list = readdir(in_dir)

Threads.@threads for i in 1:length(file_list)
    in_name = file_list[i]
    println(in_name)

    # Parse paths
    in_path = in_dir * in_name
    out_path = out_dir * replace(in_name, ".vcf" => ".parsed.tsv")
    out_path = replace(out_path, ".gz" => "")

    #@show in_name
    #@show in_path
    #@show out_path

    parse_vcfs(in_path, out_path)
end
