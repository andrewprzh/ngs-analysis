import os
import sys
import gffutils

   
def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: " + sys.argv[0] + " <GFF> <output> \n")
        exit(0)
    fn = sys.argv[1]
    db = gffutils.create_db(fn, dbfn= sys.argv[1] + '.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()
