# Define input and output file paths
input_file = "data.txt"
author_author_file = "author-author.txt"
author_conf_file = "author-conf.txt"

# Open input and output files
with open(input_file, "r") as infile, \
     open(author_author_file, "w") as author_author_out, \
     open(author_conf_file, "w") as author_conf_out:

    # Process each line in the input file
    for line in infile:
        src, dst = map(int, line.strip().split(","))  # Parse source and destination nodes
        
        # Check conditions and write to the appropriate file
        if src < 1000000000 and dst < 1000000000:
            # Both src and dst are smaller than 1000000000
            author_author_out.write(f"{src},{dst}\n")
        else:
            # Either src or dst is >= 1000000000
            if src > dst:
                src, dst = dst, src  # Ensure src is always smaller than dst
            author_conf_out.write(f"{src},{dst}\n")
