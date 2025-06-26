import argparse



def main(args: argparse.Namespace) -> None:
	...


if __name__ == "__main__":
    # create the top-level parser
    parser = argparse.ArgumentParser(description="TMAP CLI")
    parser.add_argument('--model', type=str, choices=["glide", "rwr"], required=True)
    parser.add_argument('--disease_file', type=str, required=True)
    parser.add_argument('--dataset', type=str, choices=["huri", "stringdb"], required=True)
    parser.add_argument('--cv', type=bool, default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument('--k', type=int, default=100)
    parser.add_argument('--save', type=str)

    subparsers = parser.add_subparsers(help='sub-command help')

    # RWR sub-parser
    parser_rwr = subparsers.add_parser("rwr", help='rwr help')
    parser_rwr.add_argument('--alpha', type=float, default=0.85, help='restart probability for Random Walk With Restart')

    # Glider sub-parser
    parser_glider = subparsers.add_parser('glider', help='b help')
    parser_glider.add_argument('--lamb', type=int, default=1)
    parser_glider.add_argument('--is_normalized', type=bool, default=False)
    parser_glider.add_argument('--glide_alpha', type=float, default=0.1)
    parser_glider.add_argument('--glide_beta', type=int, default=1000)
    parser_glider.add_argument('--glide_delta', type=int, default=1)
    parser_glider.add_argument('--glide_loc', type=str, choices=['cw_normalized', 'l3', 'cw'], default="cw_normalized")

args = parser.parse_args()
print(f"{args}")
main(args)