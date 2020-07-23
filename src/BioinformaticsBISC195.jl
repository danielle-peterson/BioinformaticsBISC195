module BioinformaticsBISC195

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta

# # uncomment the following line if you intend to use BioSequences types
 using BioSequences
 import BioSequences: composition, gc_content, complement, reverse_complement

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

# Your code here.
# Don't forget to export your functions!

function composition(sequence::AbstractString)
    BioSequences.composition(normalizeDNA(sequence))
end

function gc_content(sequence::AbstractString)
    BioSequences.gc_content(normalizeDNA(sequence))
end

function complement(sequence::AbstractString)
    BioSequences.complement(normalizeDNA(sequence))
end

function reverse_complement(sequence::AbstractString)
    BioSequences.reverse_complement(normalizeDNA(sequence))
end

function parse_fasta(path)
    headers = []
    sequences = []
    str = ""
    for line in eachline(path)
        if '>' ∈ line
            push!(headers, string(lstrip(line, '>')))
            if str != ""
                push!(sequences, LongDNASeq(str))
            end
            str = ""
        elseif line != ""
            for base in line
                occursin(base, "AGCTNagctn") || error("invalid base $base")
            end
            str = str * line
        end
    end
    push!(sequences, LongDNASeq(str))
    parsedfile = (headers, sequences)
    return parsedfile
end

end # module Assignment07