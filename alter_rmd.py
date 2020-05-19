import argparse



def main(infile):
    f = open(infile)
    out = open(infile.split('.')[0] + "_fixed.md", 'w')
    r_found = False
    tick_reset = False
    for line in f:
        line = line.replace("##", '')

        # Pass over the r ```r {CONTENT} ```
        if line == "```r\n" and r_found is False:
            r_found = True
            out.write(line)
        elif line == "```\n" and r_found is True:
            r_found = False
            out.write(line)

        # Reset the normal ``` {CONTENT} ```
        elif line == "```\n" and tick_reset is False:
            tick_reset = True
            out.write("<div class='r_output'>")
        elif line == "```\n" and tick_reset is True:
            tick_reset = False
            out.write("</div>")

        else:
            out.write(line)



parser = argparse.ArgumentParser(description='alter rmd.py: script to fix rmd files for better rendering',
                                 epilog='For questions or comments, please contact Matt Settles <settles@ucdavis.edu> '
                                        'or Keith Mitchell <kgmitchell@ucdavis.edu\n', add_help=True)
parser.add_argument('-i', '--input', help="Input .md rendered for Rmd to be fixed for better presentation.")

# parser.add_argument('-o', '--output', help="Output .md fixed for better presentation in template.",
#                     type=str)

options = parser.parse_args()
main(options.input)
