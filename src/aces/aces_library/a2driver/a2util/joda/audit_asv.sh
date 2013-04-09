sed -n '/\".*#.*h_ICHAR_/{			# keyword definitions
                          s/\/\*[^*]*\*\///;	# remove comments
                          s/#.*//;		# remove everything after #
                          s/.*\"//;		# remove everything up to quote
                          /[A-Z]/p;		# print keys
                         }' asv_ctl.c \
 | sort | uniq -c | sort -n
