# distribution_based_validation.py
"""
Distribution-based validation of SpatialCell AI outputs against Xenium ground truth
Compares statistical properties of gene expression across cell populations
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def validate_spatialcell_ai(validator):
    """
    Validate SpatialCell AI variants against Xenium using distribution-based metrics
    
    Parameters:
    -----------
    validator : object
        Validator object with data loading methods
    
    Returns:
    --------
    list : Results for each variant analyzed
    """
    
    print("🔄 DISTRIBUTION-BASED VALIDATION")
    print("Comparing gene expression distributions between SpatialCell AI and Xenium")
    print("=" * 60)
    
    # Load Xenium ground truth
    xenium_data = validator.load_xenium_data()
    if xenium_data is None:
        print("❌ Failed to load Xenium data")
        return None
    
    # Define SpatialCell AI variants to validate
    ai_variants = [
        ("SpatialCell AI - Visium Basic", "spatialcell_ai_visium_basic"),
        ("SpatialCell AI - Visium Advanced", "spatialcell_ai_visium_advanced"), 
        ("SpatialCell AI - HD Visium Basic 8μm", "spatialcell_ai_hdvisium_8um"),
        ("SpatialCell AI - HD Visium Basic 16μm", "spatialcell_ai_hdvisium_16um")
    ]
    
    all_results = []
    
    for variant_name, variant_path in ai_variants:
        print(f"\n{'='*60}")
        print(f"Validating: {variant_name}")
        print(f"{'='*60}")
        
        # Load SpatialCell AI variant
        ai_data = validator.load_ai_variant(variant_name, variant_path)
        if ai_data is None:
            continue
        
        # Find common genes between platforms
        ai_genes = set(ai_data.var_names)
        xenium_genes = set(xenium_data.var_names)
        common_genes = list(ai_genes & xenium_genes)
        
        print(f"Common genes: {len(common_genes)}/{len(xenium_genes)} Xenium genes")
        
        if len(common_genes) < 10:
            print(f"Insufficient common genes for {variant_name}")
            continue
        
        # Filter to common genes
        ai_filtered = ai_data[:, common_genes].copy()
        xenium_filtered = xenium_data[:, common_genes].copy()
        
        print(f"Analyzing: AI={ai_filtered.n_obs:,} cells, Xenium={xenium_filtered.n_obs:,} cells")
        
        # Calculate gene-level statistics
        gene_stats = []
        
        for gene in common_genes:
            try:
                # Extract expression values
                ai_expr = ai_filtered[:, gene].X.toarray().flatten() if hasattr(ai_filtered[:, gene].X, 'toarray') else ai_filtered[:, gene].X.flatten()
                xenium_expr = xenium_filtered[:, gene].X.toarray().flatten() if hasattr(xenium_filtered[:, gene].X, 'toarray') else xenium_filtered[:, gene].X.flatten()
                
                # Calculate statistics
                ai_mean = np.mean(ai_expr)
                xenium_mean = np.mean(xenium_expr)
                ai_std = np.std(ai_expr)
                xenium_std = np.std(xenium_expr)
                ai_detection = np.mean(ai_expr > 0)
                xenium_detection = np.mean(xenium_expr > 0)
                
                # Store statistics
                gene_stats.append({
                    'gene': gene,
                    'ai_mean': ai_mean,
                    'xenium_mean': xenium_mean,
                    'ai_std': ai_std,
                    'xenium_std': xenium_std,
                    'ai_detection_rate': ai_detection,
                    'xenium_detection_rate': xenium_detection,
                    'mean_ratio': ai_mean / xenium_mean if xenium_mean > 0 else 0,
                    'detection_ratio': ai_detection / xenium_detection if xenium_detection > 0 else 0
                })
                
            except Exception as e:
                print(f"Error processing gene {gene}: {e}")
                continue
        
        gene_stats_df = pd.DataFrame(gene_stats)
        
        if len(gene_stats_df) == 0:
            print(f"No valid genes for {variant_name}")
            continue
        
        # Calculate correlation metrics
        mean_corr = np.corrcoef(gene_stats_df['ai_mean'], gene_stats_df['xenium_mean'])[0,1]
        detection_corr = np.corrcoef(gene_stats_df['ai_detection_rate'], gene_stats_df['xenium_detection_rate'])[0,1]
        
        # Calculate ratio metrics
        mean_fold_changes = gene_stats_df['mean_ratio']
        detection_fold_changes = gene_stats_df['detection_ratio']
        
        # Filter extreme ratios for robust statistics
        valid_mean_ratios = mean_fold_changes[(mean_fold_changes > 0.01) & (mean_fold_changes < 100)]
        valid_detection_ratios = detection_fold_changes[(detection_fold_changes > 0.01) & (detection_fold_changes < 100)]
        
        median_mean_ratio = np.median(valid_mean_ratios) if len(valid_mean_ratios) > 0 else 0
        median_detection_ratio = np.median(valid_detection_ratios) if len(valid_detection_ratios) > 0 else 0
        
        # Cell-level statistics
        ai_total_counts = np.array(ai_filtered.X.sum(axis=1)).flatten()
        xenium_total_counts = np.array(xenium_filtered.X.sum(axis=1)).flatten()
        
        ai_detected_genes = np.array((ai_filtered.X > 0).sum(axis=1)).flatten()
        xenium_detected_genes = np.array((xenium_filtered.X > 0).sum(axis=1)).flatten()
        
        # Store comprehensive results
        results = {
            'variant_name': variant_name,
            'ai_cells_total': ai_data.n_obs,
            'ai_genes_total': ai_data.n_vars,
            'xenium_cells_total': xenium_data.n_obs,
            'xenium_genes_total': xenium_data.n_vars,
            'common_genes': len(common_genes),
            'genes_analyzed': len(gene_stats_df),
            
            # Gene expression correlations
            'mean_expression_correlation': mean_corr,
            'detection_pattern_correlation': detection_corr,
            
            # Distribution comparisons
            'median_mean_expression_ratio': median_mean_ratio,
            'median_detection_ratio': median_detection_ratio,
            
            # Cell-level statistics
            'ai_median_counts_per_cell': np.median(ai_total_counts),
            'xenium_median_counts_per_cell': np.median(xenium_total_counts),
            'ai_median_genes_per_cell': np.median(ai_detected_genes),
            'xenium_median_genes_per_cell': np.median(xenium_detected_genes),
            
            # Detection rates
            'ai_avg_detection': gene_stats_df['ai_detection_rate'].mean(),
            'xenium_avg_detection': gene_stats_df['xenium_detection_rate'].mean(),
            
            # Detailed statistics
            'gene_stats': gene_stats_df
        }
        
        print(f"\n📊 {variant_name} Results:")
        print(f"  Gene Expression Correlation: {mean_corr:.3f}")
        print(f"  Detection Pattern Correlation: {detection_corr:.3f}")
        print(f"  Median Expression Ratio (AI/Xenium): {median_mean_ratio:.2f}")
        print(f"  Median Detection Ratio (AI/Xenium): {median_detection_ratio:.2f}")
        print(f"  AI median counts/cell: {np.median(ai_total_counts):.1f}")
        print(f"  Xenium median counts/cell: {np.median(xenium_total_counts):.1f}")
        print(f"  Genes Successfully Analyzed: {len(gene_stats_df)}/{len(common_genes)}")
        
        all_results.append(results)
    
    # Generate validation report
    if all_results:
        generate_validation_report(all_results)
        create_validation_plots(all_results)
    
    return all_results


def generate_validation_report(all_results):
    """Generate comprehensive validation report"""
    print(f"\n{'='*80}")
    print("🎉 DISTRIBUTION-BASED VALIDATION RESULTS")
    print(f"{'='*80}")
    
    # Create comparison table
    comparison_data = []
    for result in all_results:
        comparison_data.append({
            'AI Variant': result['variant_name'].replace('SpatialCell AI - ', ''),
            'Total Cells': f"{result['ai_cells_total']:,}",
            'Common Genes': result['common_genes'],
            'Expression Corr': f"{result['mean_expression_correlation']:.3f}",
            'Detection Corr': f"{result['detection_pattern_correlation']:.3f}",
            'Expr Ratio': f"{result['median_mean_expression_ratio']:.2f}",
            'Detect Ratio': f"{result['median_detection_ratio']:.2f}",
            'AI Counts/Cell': f"{result['ai_median_counts_per_cell']:.0f}",
            'Xenium Counts/Cell': f"{result['xenium_median_counts_per_cell']:.0f}"
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    print("\n📊 VALIDATION SUMMARY TABLE:")
    print(comparison_df.to_string(index=False))
    
    # Find best performers
    best_expression = max(all_results, key=lambda x: x['mean_expression_correlation'])
    best_detection = max(all_results, key=lambda x: x['detection_pattern_correlation'])
    
    print(f"\n🏆 BEST PERFORMERS:")
    print(f"   Best Expression Correlation: {best_expression['variant_name']} (r={best_expression['mean_expression_correlation']:.3f})")
    print(f"   Best Detection Correlation: {best_detection['variant_name']} (r={best_detection['detection_pattern_correlation']:.3f})")
    
    # Save results
    comparison_df.to_csv('validation_summary.csv', index=False)
    print(f"\n💾 Results saved to: validation_summary.csv")
    
    return comparison_df


def create_validation_plots(all_results):
    """Create comprehensive validation plots"""
    print(f"\n📊 Creating validation plots...")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('SpatialCell AI Validation Against Xenium Ground Truth', 
                 fontsize=16, fontweight='bold')
    
    # Extract data for plotting
    variant_names = [r['variant_name'].replace('SpatialCell AI - ', '') for r in all_results]
    expression_corrs = [r['mean_expression_correlation'] for r in all_results]
    detection_corrs = [r['detection_pattern_correlation'] for r in all_results]
    expression_ratios = [r['median_mean_expression_ratio'] for r in all_results]
    detection_ratios = [r['median_detection_ratio'] for r in all_results]
    ai_counts = [r['ai_median_counts_per_cell'] for r in all_results]
    xenium_counts = [r['xenium_median_counts_per_cell'] for r in all_results]
    
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4']
    
    # Plot 1: Expression correlations
    bars1 = axes[0,0].bar(variant_names, expression_corrs, color=colors[:len(variant_names)], alpha=0.8)
    axes[0,0].set_ylabel('Correlation Coefficient')
    axes[0,0].set_title('Gene Expression Correlation')
    axes[0,0].set_ylim(0, 1)
    axes[0,0].axhline(y=0.7, color='red', linestyle='--', alpha=0.5, label='r=0.7')
    
    for bar, corr in zip(bars1, expression_corrs):
        axes[0,0].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                      f'{corr:.3f}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 2: Detection correlations
    bars2 = axes[0,1].bar(variant_names, detection_corrs, color=colors[:len(variant_names)], alpha=0.8)
    axes[0,1].set_ylabel('Correlation Coefficient')
    axes[0,1].set_title('Detection Pattern Correlation')
    axes[0,1].set_ylim(0, 1)
    axes[0,1].axhline(y=0.8, color='red', linestyle='--', alpha=0.5, label='r=0.8')
    
    for bar, corr in zip(bars2, detection_corrs):
        axes[0,1].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                      f'{corr:.3f}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 3: Expression ratios
    bars3 = axes[0,2].bar(variant_names, expression_ratios, color=colors[:len(variant_names)], alpha=0.8)
    axes[0,2].set_ylabel('Ratio (AI/Xenium)')
    axes[0,2].set_title('Median Expression Level Ratio')
    axes[0,2].axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Perfect match')
    
    for bar, ratio in zip(bars3, expression_ratios):
        axes[0,2].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.05,
                      f'{ratio:.2f}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 4: Detection ratios
    bars4 = axes[1,0].bar(variant_names, detection_ratios, color=colors[:len(variant_names)], alpha=0.8)
    axes[1,0].set_ylabel('Ratio (AI/Xenium)')
    axes[1,0].set_title('Median Detection Rate Ratio')
    axes[1,0].axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Perfect match')
    
    for bar, ratio in zip(bars4, detection_ratios):
        axes[1,0].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.05,
                      f'{ratio:.2f}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 5: Counts per cell comparison
    x_pos = np.arange(len(variant_names))
    width = 0.35
    
    bars5a = axes[1,1].bar(x_pos - width/2, ai_counts, width, label='SpatialCell AI', color='#3498db', alpha=0.8)
    bars5b = axes[1,1].bar(x_pos + width/2, xenium_counts, width, label='Xenium', color='#e74c3c', alpha=0.8)
    
    axes[1,1].set_ylabel('Median Counts per Cell')
    axes[1,1].set_title('Counts per Cell Comparison')
    axes[1,1].set_xticks(x_pos)
    axes[1,1].set_xticklabels(variant_names)
    axes[1,1].legend()
    
    # Plot 6: Composite scores
    composite_scores = []
    for i in range(len(all_results)):
        # Normalize correlations (higher is better)
        expr_score = expression_corrs[i] if not np.isnan(expression_corrs[i]) else 0
        detect_score = detection_corrs[i] if not np.isnan(detection_corrs[i]) else 0
        
        # Normalize ratios (closer to 1 is better)
        expr_ratio_score = 1 - abs(expression_ratios[i] - 1) if expression_ratios[i] > 0 else 0
        detect_ratio_score = 1 - abs(detection_ratios[i] - 1) if detection_ratios[i] > 0 else 0
        
        # Ensure scores are between 0 and 1
        expr_ratio_score = max(0, min(1, expr_ratio_score))
        detect_ratio_score = max(0, min(1, detect_ratio_score))
        
        # Composite score (geometric mean)
        composite = (expr_score * detect_score * expr_ratio_score * detect_ratio_score) ** 0.25
        composite_scores.append(composite)
    
    bars6 = axes[1,2].bar(variant_names, composite_scores, color=colors[:len(variant_names)], alpha=0.8)
    axes[1,2].set_ylabel('Composite Score')
    axes[1,2].set_title('Overall Performance Score')
    axes[1,2].set_ylim(0, 1)
    
    for bar, score in zip(bars6, composite_scores):
        axes[1,2].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                      f'{score:.3f}', ha='center', va='bottom', fontweight='bold')
    
    # Rotate x-labels for better readability
    for ax in axes.flat:
        ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('validation_results.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("📊 Validation plots saved as 'validation_results.png'")
    
    return fig


# Main execution
if __name__ == "__main__":
    print("🚀 Running SpatialCell AI Validation")
    print("Distribution-based comparison with Xenium ground truth")
    
    # Initialize validator (implement based on your data structure)
    # validator = SpatialCellValidator(data_path="path/to/data")
    # results = validate_spatialcell_ai(validator)
