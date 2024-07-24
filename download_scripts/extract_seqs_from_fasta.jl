using FASTX
using ThreadSafeDicts

function extract_seq(filename::AbstractString, headers_path::AbstractString)
    headers = readlines(headers_path)
    
    all_headers = Dict{String, Int64}()
    fasta_index = 1

    FASTA.Reader(open(filename)) do reader
        for record in reader
            header = FASTA.identifier(record)
            all_headers[header] = fasta_index
            fasta_index += 1
            println(fasta_index)
        end
    end

    return all_headers
end


# Example usage
fasta_file = "/mnt/c/git_repos/wuhu_rooting/data/genomes/gisaid/gisaid_sequences.230324.fna"
headers_file = "/mnt/c/git_repos/wuhu_rooting/data/metadata/gisaid_metatadata.220324.1Apr20.accessions_only.txt"
out_filt = "test.fna"

extract_seq(fasta_file, headers_file)
