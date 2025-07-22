def merge_supervised_cluster_rankings(rankings, disease_genes):
        max_score = max(ranking[0][1] for ranking in rankings)

        def assign_label(score, threshold):
                return score / max_score if score >= threshold else score

        def get_threshold(ranking, min_n=20, fallback_q=75, visualize=True):
                # valley-finding heuristic

                # Parameters
                # ----------
                # ranking : list[(gene, score)]
                # min_n   : int        # require at least this many scores to attempt valley search
                # fallback_q : float   # percentile to return when valley search fails
                
                # --- 0. collect scores and exit early if we have nothing ---
                scores = np.asarray([s for gene, s in ranking if 0 < s <= 1.0 and gene not in disease_genes], dtype=float)
                scores = scores[np.isfinite(scores)]
                if scores.size == 0:
                        return 0

                # --- 1. skip valley-finding on tiny clusters; use quantile fallback ---
                if scores.size < min_n:
                        return float(np.percentile(scores, fallback_q))

                # --- 2. histogram & smoothing (unchanged except for adaptive sigma) ---
                counts, bin_edges = np.histogram(scores, bins=20)
                log_counts = np.where(counts > 0, np.log10(counts), np.nan)

                sigma = max(1, log_counts.size // 50)           # light adaptive smoothing
                smoothed = gaussian_filter1d(log_counts, sigma=sigma)

                minima = argrelextrema(smoothed, np.less)[0]
                maxima = argrelextrema(smoothed, np.greater)[0]

                if visualize:
                        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                        plt.figure(figsize=(8, 5))
                        plt.bar(bin_centers, log_counts, width=np.diff(bin_edges), alpha=0.4, label="Histogram")
                        plt.plot(bin_centers, smoothed, label="Smoothed")
                        plt.scatter(bin_centers[minima], smoothed[minima], color='red', label='Minima', zorder=3)
                        plt.scatter(bin_centers[maxima], smoothed[maxima], color='green', label='Maxima', zorder=3)
                        plt.xlabel("Score")
                        plt.ylabel("Frequency (log scale)")
                        plt.title("Histogram + Smoothed Curve (for valley detection)")
                        plt.legend()
                        plt.show()

                # --- 3. find valley with widest *score-space* span ---
                widest_span = 0.0
                best_min = None
                for m in minima:
                        left = maxima[maxima < m]
                        right = maxima[maxima > m]
                        if left.size == 0 or right.size == 0:
                                continue

                        span = bin_edges[right.min() + 1] - bin_edges[left.max()]
                        if span > widest_span:
                                widest_span = span
                                best_min = m

                # --- 4. return threshold or fallback ---
                if best_min is not None:
                        return 0.5 * (bin_edges[best_min] + bin_edges[best_min + 1])

                return float(np.percentile(scores, fallback_q))

        
        threshold = sum(get_threshold(ranking) for ranking in rankings) / len(rankings)
        print("Threshold:", threshold)

        final_scores = {}
        for ranking in rankings:
                for gene, score in ranking:
                        final_scores[gene] = final_scores.get(gene, 0) + assign_label(score, threshold)  # find general or adaptive threshold for other diseases

        return sorted(list(final_scores.items()), key=lambda x: -x[1])
        """
        max_score = max(ranking[0][1] for ranking in rankings)

        for ranking in rankings:
                for i in range(len(ranking)):
                        ranking[i] = (ranking[i][0], pow(ranking[i][1] / max_score, 4))

        from collections import defaultdict
        max_rankings = defaultdict(int)

        for ranking in rankings:
               for gene, score in ranking:
                        max_rankings[gene] += score 
        
        return sorted(list(max_rankings.items()), key=lambda x: -x[1])
        """

def lenore_merging(og_ranking, rankings):
        """
        go through og ranking
        for each place, check if it appears at that ranking or higher in any cluster ranking
        if not, lower the score by some number
        re-sort og ranking
        """
        best_placements = defaultdict(int)
        for ranking in rankings:
                for i, (gene, score) in enumerate(ranking):
                        best_placements[gene] = min(best_placements[gene], i)

        final_scores = {}
        for i, (gene, score) in enumerate(og_ranking):
                if best_placements[gene] <= i:
                        final_scores[gene] = score * 1.2

                else:
                        final_scores[gene] = score

        return sorted(list(final_scores.items()), key=lambda x: -x[1])


def main():
        rankings_path = input("Enter path to rankings folder: ")
        rankings


if __name__ == "__main__":
        args = parser.parse_args()
        main(args.network, args.genelist, args.out)