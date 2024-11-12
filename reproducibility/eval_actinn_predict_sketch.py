import pandas as pd
import argparse

def eval(path):
    pred_df = pd.read_csv(path, sep='\t', header=0, index_col=0)
    pred_df.columns = ['pred_label']
    true_df = pd.read_csv('actinn_data/test.age.label.txt', sep='\t', header=0, index_col=0)
    true_df.columns = ['true_label']
    df = true_df.join(pred_df)
    acc = sum(df['pred_label']==df['true_label']) / df.shape[0]
    print('Test ACC=', acc)
    return acc

def get_parser(parser=None):
    if parser == None:
        parser = argparse.ArgumentParser()
    parser.add_argument("-sm", "--method", type=str, help="Sketch method.")
    parser.add_argument("-ss", "--sketch_size", type=int, help="Sketch_size.")
    parser.add_argument("-rs", "--random_seed", type=int, help="Random seed.")
    return parser

if __name__ == '__main__':

    parser = get_parser()
    args = parser.parse_args()

    path = "experiments/seed%d/predicted_label.%s.%d.txt" % (args.random_seed, args.method, args.sketch_size)
    eval(path)
