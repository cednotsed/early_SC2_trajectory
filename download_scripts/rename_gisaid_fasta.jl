using FASTX
using ThreadSafeDicts

# Function to read FASTA file and extract sequences based on headers
function parse_headers(fasta_file::AbstractString, out_file::AbstractString)
    reader = open(FASTA.Reader, fasta_file)
    
    # Open FASTA file and iterate over records using multiple threads
    counter = 0
    io = open(out_file, "a")

    for record in reader
        new_header = split(FASTA.identifier(record), '|')[1]
        counter += 1
        println(io, ">$new_header")
        println(io, String(FASTA.sequence(record)))
        println(counter)
    end
    
    close(reader)
    close(io)
end

# Example usage
fasta_file = "/mnt/c/git_repos/wuhu_rooting/data/genomes/gisaid/gisaid_sequences.230324.fna"
out_file = "/mnt/c/git_repos/wuhu_rooting/data/genomes/gisaid/gisaid_sequences.230324.renamed.fna"

parse_headers(fasta_file, out_file)
