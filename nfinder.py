### Finds and reports regions of N's in FASTA records
# Coded for python 3.7
import sys, re, os.path, argparse, io

#### Help string:
# usage = """%s: Finds and reports regions of N's in FASTA records.
# JI Ohlsson, UTK 2020, johlsson@utk.edu
#
# Usage: %s <FASTA file> [ options ]
#
# Options:
#    -h		show this help
# """ % (sys.argv[0], sys.argv[0])

### Get arguments
parser = argparse.ArgumentParser(
	description="""Finds and reports regions of N's in FASTA records.
I Ohlsson, UTK 2020, johlsson@utk.edu""")
parser.add_argument('filename', metavar='<FASTA file>', help='FASTA sequence file')
parser.add_argument('-s', '--split', type=int, metavar='X', nargs='?', const=10, default=None, help='Split sequences around N region of length X (default=10)')
parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), metavar='F', nargs=1, default=None, help='Output split sequences to file F')
args = parser.parse_args()

# If -s is defined without -o, will append split seqs after output
if args.outfile is None:
	outstream = io.StringIO()
else:
	outstream = args.outfile[0]
#sys.exit(args)

### Prepare regular expressions
# NOTE: remember to uppercase string before matching
# NOTE: considers as "sequence" only characters in the official
# 	FASTA nucleic acid codes
hed_re = re.compile("^>")	# Matches FASTA header
seq_re = re.compile("^[ACGTURYKMSWBDHVN-]+")	# Matches FASTA sequence
#seq_strict_re = re.compile("^[ACGTUN-]+")	# Matches only non-ambiguous NA, plus N(any) and -(gap)
n_re = re.compile("(N+)")

def report_ns(pattern, header, sequence):
	"""Search for N regions in sequence and report them.
	header: nlength nstart nend (for each match)
	"""
	# Get an iterator containing all (non-overlapping) matches for pattern:
	n_regions = list(pattern.finditer(sequence))
	
	if n_regions:
		# List all N regions
		for r in n_regions:
			(start, end) = r.span(0)
			print("\t".join([header, str(end - start), str(start), str(end)]))
	else:
		# No N regions found!
		print("\t".join([header, "-", "-", "-"]))

def split_seq(length, header, sequence, stream: io.TextIOBase):
	"""Split sequence around N regions of given minimum length.
	Sub-sequence headers are appended with "_part_#_<start-end>"
	NOTE: assumes there is an outfile to write to!
	"""
	# Generate a pattern:
	pattern_string = "N{%d,}" % length
	pattern = re.compile(pattern_string)
	
	# Get an iterator containing all (non-overlapping) matches for pattern:
	n_regions = list(pattern.finditer(sequence))
	
	last_end = 0	# End of last match; starts at beginning of sequence
	count = 0
	if n_regions:
		# Split sequence around N regions of specified length
		for r in n_regions:
			# Get coordinates of N region
			(start, end) = r.span(0)
			# Get subsequence from end of last match to start of this
			subseq = sequence[last_end:start]
			# Avoid printing empty sequences:
			if len(subseq) > 0:
				# Print subsequence FASTA record
				# NOTE: uses a TextIOWrapper, provided by argparser, to write out
				# NOTE: ASSUMES THERE IS A FILE TO WRITE TO
				stream.write("%s_part_%d_coords_%d-%d_newlen_%d\n" % (header, count, last_end, start, start-last_end))
				stream.write(subseq + "\n")
				# Update subsequence count
				count += 1
			last_end = end
			
		# Write final subsequence (end of last match to end of sequence)
		subseq = sequence[last_end:]
		# Again, avoid printing empty subsequences:
		if len(subseq) > 0:
			stream.write("%s_part_%d_coords_%d-%d_newlen_%d\n" % (header, count, last_end, len(sequence), len(sequence)-last_end))
			stream.write(subseq + "\n")
	else:
		# No sufficiently long N regions: just print sequence
		stream.write(header + "\n")
		stream.write(sequence + "\n")
		pass
	
	

### Process FASTA file
# Make sure each sequence is on a single line, to allow safe
# searching for N's across line breaks
header = None
sequence = ""
linenum = 0

with open(args.filename) as file:
	print("Record\tLength\tStart\tEnd")
	for line in file:
		line = line.rstrip()
		if line.isspace():
			# Skip empty lines
			continue
		if hed_re.match(line):
			# If line is a header:
			# Process previous record
			if header is not None:
				report_ns(n_re, header, sequence)
				if args.split:
					# Split sequence if asked to:
					split_seq(args.split, header, sequence, outstream)
			# Save new header and reset sequence
			header = line
			sequence = ""
		elif seq_re.fullmatch(line.upper()):
			# If line is sequence:
			# TODO: option to use strict sequence matching?
			# Add to current sequence
			sequence += line
		else:
			# Couldn't match normal FASTA lines! Stop!
			# We don't want to continue with parts left out of a sequence.
			sys.exit("""ERROR: line %d in %s does not look like FASTA:\n
						%s""" % (linenum, args.filename, line))
		
		# Update line number
		linenum += 1
	
	# After iteration, process final record
	if header is not None:
		report_ns(n_re, header, sequence)
		if args.split:
			# Split sequence if asked to:
			split_seq(args.split, header, sequence, outstream)

# Finally, if not writing to a file, dump split sequences to terminal:
if args.split and not args.outfile:
	if os.path.getsize(args.filename) > 4096:
		# Try to discourage dumping large files to stdout:
		sys.stderr.write("WARNING: You did not specify an output file (-o); split sequences will be output to the terminal.\nContinue? (y/n)\n")
		if not input().lower().strip()[:1] == "y":
			sys.exit(0)
	# Print to terminal
	outstream.seek(0)
	for line in outstream:
		print(line, end='')