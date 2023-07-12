import argparse
import amr_functions as amr
import os

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Python Program")

    parser.add_argument("-p", "--patheda", default='lib\\genes_eda_data',
                        help="The EDA path argument (default: 'lib\\genes_eda_data' in the current directory)")
    parser.add_argument("-pr", "--pathres", default='results',
                        help="The path argument (default: 'results' in the current directory))")
    parser.add_argument("-k", "--kmer", type=int, choices=range(2, 12), default=8,
                        help="The k_mer argument (an integer from 2 to 11, default: 5)")
    
    args = parser.parse_args()
    
    print("The value of the path eda argument is:", args.patheda)
    print("The value of the path results argument is:", args.pathres)
    print("The value of the kmer argument is:", args.kmer)


    titles = {}
    titles = amr.read_ascii_art(os.path.join(os.getcwd(), 'titles.txt'))
    print(titles['results'])

    # Process the folders in the specified base path and generate dataframes and CSV files.
    amr.compile_results(args.pathres,args.patheda,args.kmer)
    # Read the CSV file, extract column names, and plot scatter points.
    amr.plot_scatter_full(args.patheda)
    # Analyze the data from a CSV file and generate visualizations and summary tables.
    amr.plot_scatter_top(args.patheda)



if __name__ == "__main__":
    main()
