import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
import pingouin as pg
from matplotlib.patches import Patch
from scipy.stats import shapiro, levene
import matplotlib.ticker as mtick

# Read the CSV file
df = pd.read_csv('cfu count thesis.csv', index_col=0)

# Data preprocessing
# Extract the position (T/M/B) and treatment (con/bot/pb) from the index
data = []

for idx in df.index:
    position = idx[0]  # First character (T, M, or B)
    treatment = idx[1:].replace('1', '').replace('2', '').replace('3', '')  # Extract treatment
    rep_number = int(idx[-1])  # Repetition number
    
    for col in range(1, 6):  # Columns 1-5 contain the actual measurements
        value = df.iloc[df.index.get_loc(idx), col-1]
        data.append({
            'Position': position,
            'Treatment': treatment,
            'Replicate': rep_number,
            'Measurement': col,
            'CFU': value
        })

# Create a tidy dataframe
tidy_df = pd.DataFrame(data)

# Map positions and treatments to more readable labels
position_map = {'T': 'Top', 'M': 'Middle', 'B': 'Bottom'}
treatment_map = {'con': 'Control', 'bot': 'Botector', 'pb': 'Potassium Bicarbonate'}
tidy_df['Position'] = tidy_df['Position'].map(position_map)
tidy_df['Treatment'] = tidy_df['Treatment'].map(treatment_map)

# Create merged treatment-position column for some visualizations
tidy_df['Treatment_Position'] = tidy_df['Treatment'] + '_' + tidy_df['Position']

# Calculate summary statistics
summary_stats = tidy_df.groupby(['Treatment', 'Position'])['CFU'].agg(['count', 'mean', 'std', 'min', 'max', 'sem']).reset_index()
summary_stats['cv'] = (summary_stats['std'] / summary_stats['mean']) * 100  # coefficient of variation

# Calculate relative efficacy compared to control
control_means = summary_stats[summary_stats['Treatment'] == 'Control'].set_index('Position')['mean']
for idx, row in summary_stats.iterrows():
    if row['Treatment'] != 'Control':
        control_mean = control_means[row['Position']]
        summary_stats.at[idx, 'efficacy'] = ((control_mean - row['mean']) / control_mean) * 100

# Normality test
normality_results = pd.DataFrame(columns=['Treatment_Position', 'W', 'p', 'Normal'])
for group in tidy_df['Treatment_Position'].unique():
    subset = tidy_df[tidy_df['Treatment_Position'] == group]['CFU']
    stat, p = shapiro(subset)
    normality_results = pd.concat([normality_results, pd.DataFrame({
        'Treatment_Position': [group], 
        'W': [stat], 
        'p': [p], 
        'Normal': [p > 0.05]
    })])

# Homogeneity of variance test
levene_stat, levene_p = levene(*[group['CFU'].values for name, group in tidy_df.groupby('Treatment_Position')])
equal_variance = levene_p > 0.05

# Two-way ANOVA
model = ols('CFU ~ C(Treatment) + C(Position) + C(Treatment):C(Position)', data=tidy_df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

# Post-hoc tests
# For Treatment
posthoc_treatment = pg.pairwise_tukey(data=tidy_df, dv='CFU', between='Treatment')

# For Position
posthoc_position = pg.pairwise_tukey(data=tidy_df, dv='CFU', between='Position')

# For interaction
posthoc_interaction = []
for pos in tidy_df['Position'].unique():
    subset = tidy_df[tidy_df['Position'] == pos]
    if len(subset) > 0:
        result = pg.pairwise_tukey(data=subset, dv='CFU', between='Treatment')
        result['Position'] = pos
        posthoc_interaction.append(result)
posthoc_interaction_df = pd.concat(posthoc_interaction)

# Create visualizations
plt.style.use('ggplot')
colors = {'Control': '#e41a1c', 'Botector': '#377eb8', 'Potassium Bicarbonate': '#4daf4a'}

# 1. Enhanced box plot
plt.figure(figsize=(14, 8))
ax = sns.boxplot(x='Position', y='CFU', hue='Treatment', data=tidy_df, palette=colors)

# Add individual data points
sns.stripplot(x='Position', y='CFU', hue='Treatment', data=tidy_df, 
              size=4, alpha=0.6, dodge=True, palette=colors, jitter=True)

plt.title('CFU Counts by Position and Treatment', fontsize=16, fontweight='bold')
plt.ylabel('Colony Forming Units (CFU)', fontsize=14)
plt.xlabel('Position on Slope', fontsize=14)
plt.legend(title='Treatment', title_fontsize=12, fontsize=10, loc='upper right')
plt.grid(True, linestyle='--', alpha=0.7)

# Add statistical significance annotations
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:3], labels[:3], title='Treatment')

plt.tight_layout()
plt.savefig('enhanced_cfu_boxplot.png', dpi=300)
plt.close()

# 2. Enhanced bar plot with error bars
fig, ax = plt.subplots(figsize=(14, 8))

# Prepare data for grouped bar plot
summary = tidy_df.groupby(['Position', 'Treatment'])['CFU'].agg(['mean', 'sem']).reset_index()

# Define the x locations for the groups
positions = np.arange(len(position_map))
bar_width = 0.25
opacity = 0.8

# Plot each treatment as grouped bars
for i, treatment in enumerate(['Control', 'Botector', 'Potassium Bicarbonate']):
    treatment_data = summary[summary['Treatment'] == treatment]
    means = [treatment_data[treatment_data['Position'] == pos]['mean'].values[0] if not treatment_data[treatment_data['Position'] == pos].empty else 0 for pos in position_map.values()]
    errors = [treatment_data[treatment_data['Position'] == pos]['sem'].values[0] if not treatment_data[treatment_data['Position'] == pos].empty else 0 for pos in position_map.values()]
    
    ax.bar(positions + (i-1)*bar_width, means, bar_width,
           alpha=opacity, color=colors[treatment], label=treatment,
           yerr=errors, capsize=5)

# Add labels, title and legend
ax.set_xlabel('Position on Slope', fontsize=14)
ax.set_ylabel('Mean CFU Count', fontsize=14)
ax.set_title('Mean CFU Counts by Position and Treatment', fontsize=16, fontweight='bold')
ax.set_xticks(positions)
ax.set_xticklabels(position_map.values())
ax.legend(title='Treatment')

plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('grouped_cfu_barplot.png', dpi=300)
plt.close()

# 3. Heat map of CFU counts
plt.figure(figsize=(12, 8))
heatmap_data = tidy_df.groupby(['Position', 'Treatment'])['CFU'].mean().unstack()
sns.heatmap(heatmap_data, annot=True, fmt=".4f", cmap="YlOrRd", linewidths=.5)
plt.title('Mean CFU Counts Heatmap', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('cfu_heatmap.png', dpi=300)
plt.close()

# 4. Violin plot with embedded boxplot
plt.figure(figsize=(14, 8))
ax = sns.violinplot(x='Position', y='CFU', hue='Treatment', data=tidy_df, 
                    palette=colors, split=False, inner='box', 
                    linewidth=1)
plt.title('Distribution of CFU Counts by Position and Treatment', fontsize=16, fontweight='bold')
plt.ylabel('Colony Forming Units (CFU)', fontsize=14)
plt.xlabel('Position on Slope', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('cfu_violin_plot.png', dpi=300)
plt.close()

# 5. Efficacy comparison chart
efficacy_data = summary_stats[summary_stats['Treatment'] != 'Control'].copy()
plt.figure(figsize=(12, 7))
bars = sns.barplot(x='Position', y='efficacy', hue='Treatment', data=efficacy_data, 
                 palette={'Botector': colors['Botector'], 'Potassium Bicarbonate': colors['Potassium Bicarbonate']})

# Add percentage signs to y-axis
ax = plt.gca()
ax.yaxis.set_major_formatter(mtick.PercentFormatter())

# Add value labels on the bars
for container in bars.containers:
    bars.bar_label(container, fmt='%.4f%%')

plt.title('Treatment Efficacy Compared to Control', fontsize=16, fontweight='bold')
plt.ylabel('Reduction in CFU (%)', fontsize=14)
plt.xlabel('Position on Slope', fontsize=14)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('treatment_efficacy.png', dpi=300)
plt.close()

# 6. Interaction plot
plt.figure(figsize=(12, 7))
interaction = tidy_df.groupby(['Position', 'Treatment'])['CFU'].mean().reset_index()
interaction_pivot = interaction.pivot(index='Position', columns='Treatment', values='CFU')

for treatment, color in colors.items():
    plt.plot(interaction_pivot.index, interaction_pivot[treatment], marker='o', 
             linewidth=2, markersize=8, label=treatment, color=color)

plt.title('Interaction Plot: Treatment × Position', fontsize=16, fontweight='bold')
plt.ylabel('Mean CFU Count', fontsize=14)
plt.xlabel('Position on Slope', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(title='Treatment')
plt.tight_layout()
plt.savefig('enhanced_interaction_plot.png', dpi=300)
plt.close()

# Generate markdown report
markdown = f"""# Colony-Forming Unit (CFU) Analysis Report

## Overview

This analysis examines the effects of three treatments (Control, Botector, and Potassium Bicarbonate) across three positions on a slope (Top, Middle, Bottom) on Colony-Forming Unit (CFU) counts. The goal is to determine if there are significant differences between treatments and if position on the slope affects CFU counts.

### Experimental Context

Colony-Forming Units (CFUs) are a measure of viable fungal or bacterial cells in a sample. In this experiment, three treatments were applied across three positions on a slope:

- **Treatments**: Control (no treatment), Botector (a biological control agent), and Potassium Bicarbonate (a chemical fungicide)
- **Positions**: Top, Middle, and Bottom of slope
- **Replication**: Each treatment-position combination was replicated three times, with five measurements taken per replicate

## Summary Statistics

The table below presents key summary statistics for each combination of treatment and position:

```
{summary_stats.to_markdown(index=False, floatfmt=".4f")}
```

*Note: 'count' represents number of measurements, 'mean' and 'std' are average and standard deviation of CFU counts, 'sem' is standard error of mean, 'cv' is coefficient of variation (%), and 'efficacy' shows percent reduction compared to control.*

## Statistical Analysis

Before conducting the main analysis, we need to check if our data meets the necessary assumptions for parametric testing.

### Assumptions Testing

#### Normality Test (Shapiro-Wilk)

The Shapiro-Wilk test examines whether data follows a normal distribution. This is important because Analysis of Variance (ANOVA) assumes that the data within each group is normally distributed.

- **W statistic**: Values closer to 1 indicate normality
- **p-value**: If p > 0.05, we cannot reject the assumption of normality
- **Interpretation**: Data is considered normally distributed when p > 0.05

```
{normality_results.to_markdown(index=False, floatfmt=".4f")}
```

#### Homogeneity of Variance (Levene's Test)

Levene's test examines whether the different groups have similar variances. ANOVA assumes that all groups have similar spread of data.

- **Test statistic**: {levene_stat:.4f}
- **p-value**: {levene_p:.4f}
- **Equal variance assumption**: {'Met' if equal_variance else 'Not met'}
- **Interpretation**: {'Since p > 0.05, the groups have similar variances, meeting this ANOVA assumption.' if equal_variance else 'Since p < 0.05, the groups have different variances, violating this ANOVA assumption. Results should be interpreted with caution.'}

### Two-Way ANOVA Results

**What is ANOVA?** Analysis of Variance (ANOVA) is a statistical method used to compare means of multiple groups. In this case, we're using a two-way ANOVA because we have two factors: Treatment and Position. This analysis helps us understand:

1. **Main effect of Treatment**: Does the type of treatment (Control, Botector, or Potassium Bicarbonate) significantly affect CFU counts, regardless of position?
2. **Main effect of Position**: Does the position on the slope (Top, Middle, or Bottom) significantly affect CFU counts, regardless of treatment?
3. **Interaction effect**: Do treatments behave differently depending on position? For example, does Botector work better at the Top position compared to other positions?

**How to interpret the results**:
- The **sum_sq** column shows the sum of squares, or variation explained by each factor
- The **F** value is the test statistic - higher values indicate a stronger effect
- The **PR(>F)** column shows the p-value - values less than 0.05 indicate statistical significance

```
{anova_table.round(4).to_markdown()}
```

#### Key Findings from ANOVA:
- **Treatment effect**: {'Significant' if anova_table.loc['C(Treatment)', 'PR(>F)'] < 0.05 else 'Not significant'} (p = {anova_table.loc['C(Treatment)', 'PR(>F)']:.4f}) - {'This means that the different treatments do have significantly different effects on CFU counts.' if anova_table.loc['C(Treatment)', 'PR(>F)'] < 0.05 else 'This means that we cannot conclude that the treatments have different effects on CFU counts.'}
- **Position effect**: {'Significant' if anova_table.loc['C(Position)', 'PR(>F)'] < 0.05 else 'Not significant'} (p = {anova_table.loc['C(Position)', 'PR(>F)']:.4f}) - {'This means that the position on the slope significantly affects CFU counts.' if anova_table.loc['C(Position)', 'PR(>F)'] < 0.05 else 'This means that the position on the slope does not significantly affect CFU counts.'}
- **Interaction effect**: {'Significant' if anova_table.loc['C(Treatment):C(Position)', 'PR(>F)'] < 0.05 else 'Not significant'} (p = {anova_table.loc['C(Treatment):C(Position)', 'PR(>F)']:.4f}) - {'This means that the effect of treatments depends on the position on the slope. In other words, some treatments may work better at certain positions than others.' if anova_table.loc['C(Treatment):C(Position)', 'PR(>F)'] < 0.05 else 'This means that treatments have a consistent effect regardless of position on the slope.'}

### Post-hoc Tests

When ANOVA indicates significant differences, we need to perform follow-up tests to determine exactly which groups differ from each other. This is where Tukey's Honest Significant Difference (HSD) test comes in.

#### Tukey's HSD for Treatment Factor

**What is Tukey's HSD?** This test compares all possible pairs of treatments to identify which specific treatments differ significantly from each other. It controls for the increased risk of false positives when making multiple comparisons.

**How to interpret the results**:
- **A** and **B** columns show which treatments are being compared
- **mean(A)** and **mean(B)** show the average CFU counts for each treatment
- **diff** shows the difference between these means
- **p-tukey** is the adjusted p-value - values less than 0.05 indicate a significant difference
- **hedges** is an effect size measure

```
{posthoc_treatment.to_markdown(index=False, floatfmt=".4f")}
```

#### Tukey's HSD for Position Factor

Similar to the treatment comparison, this test identifies which specific positions (Top, Middle, Bottom) differ significantly from each other in terms of CFU counts.

```
{posthoc_position.to_markdown(index=False, floatfmt=".4f")}
```

#### Treatment Comparison Within Each Position

Since we're also interested in how treatments perform at each specific position, these tests compare treatments separately within each position group. This helps us understand if, for example, Botector is particularly effective at the Top position but not at others.

```
{posthoc_interaction_df.to_markdown(index=False, floatfmt=".4f")}
```

## Visualizations

Each visualization presents a different perspective on the data to help understand the relationships between treatments and positions.

### 1. Enhanced Box Plot with Individual Data Points
![Enhanced Box Plot](enhanced_cfu_boxplot.png)

**What this shows:** This plot displays the distribution of CFU counts for each treatment across the three positions on the slope. 

**How to interpret:**
- The colored boxes show the interquartile range (middle 50% of data) for each treatment
- The horizontal line inside each box represents the median value
- The "whiskers" extend to the minimum and maximum values (excluding outliers)
- Individual points represent actual measurements, allowing you to see the raw data
- The x-axis groups data by position (Top, Middle, Bottom), with different colors representing the three treatments

**Key insights:** This visualization helps identify differences in both the central tendency and the spread of CFU counts across treatments and positions, while showing the actual data points.

### 2. Grouped Bar Plot
![Grouped Bar Plot](grouped_cfu_barplot.png)

**What this shows:** This plot presents the mean CFU counts for each treatment at each position, with error bars indicating standard error of the mean.

**How to interpret:**
- Each position (Top, Middle, Bottom) has three bars representing the three treatments
- The height of each bar represents the mean CFU count
- Error bars show the standard error, giving an indication of the reliability of the mean
- Different colors distinguish between the three treatments

**Key insights:** This visualization clearly shows the average performance of each treatment at each position, making it easy to compare treatment efficacy across positions.

### 3. Heat Map of CFU Counts
![Heat Map](cfu_heatmap.png)

**What this shows:** This heat map uses color intensity to visualize mean CFU counts for each treatment-position combination.

**How to interpret:**
- Rows represent positions on the slope
- Columns represent treatments
- Color intensity corresponds to CFU count (darker red = higher count)
- The numbers in each cell show the exact mean CFU count

**Key insights:** Heat maps provide a quick visual overview of which treatment-position combinations have the highest and lowest CFU counts, making patterns easier to spot.

### 4. Violin Plot
![Violin Plot](cfu_violin_plot.png)

**What this shows:** Violin plots combine aspects of box plots with density plots to show the distribution of CFU counts.

**How to interpret:**
- The width of each "violin" at any point represents the density of data at that CFU count
- Wider sections indicate more data points at that CFU value
- Inside each violin is a small box plot showing median and interquartile range
- The x-axis groups by position, with colors distinguishing treatments

**Key insights:** This visualization reveals the full distribution shape of each dataset, showing whether CFU counts are concentrated in certain ranges or more evenly distributed.

### 5. Treatment Efficacy Compared to Control
![Treatment Efficacy](treatment_efficacy.png)

**What this shows:** This bar chart displays the percentage reduction in CFU counts for each treatment compared to the control at each position.

**How to interpret:**
- The x-axis shows the three positions
- The y-axis shows percent reduction in CFU counts compared to the control
- Higher percentages indicate greater efficacy (more reduction in CFUs)
- Different colors distinguish between treatments
- The percentage values are labeled on each bar

**Key insights:** This visualization directly shows how effective each treatment is at reducing CFU counts compared to the control at each position on the slope.

### 6. Enhanced Interaction Plot
![Enhanced Interaction Plot](enhanced_interaction_plot.png)

**What this shows:** This interaction plot illustrates how the effect of treatments varies across positions.

**How to interpret:**
- The x-axis represents position on the slope
- The y-axis shows mean CFU counts
- Each line represents a different treatment
- Non-parallel lines indicate an interaction effect (the treatment effect varies by position)

**Key insights:** This visualization helps identify if treatments perform consistently across positions or if their effectiveness depends on position. When lines cross or have very different slopes, it suggests that position affects how well a treatment works.

## Conclusion

""" 

# Add conclusion based on actual results
if anova_table.loc['C(Treatment)', 'PR(>F)'] < 0.05:
    significant_treatments = posthoc_treatment[posthoc_treatment['p-tukey'] < 0.05]
    if not significant_treatments.empty:
        markdown += "The analysis reveals a significant effect of treatment on CFU counts. "
        for _, row in significant_treatments.iterrows():
            markdown += f"The {row['A']} and {row['B']} treatments differ significantly (p = {row['p-tukey']:.4f}). "
    else:
        markdown += "Although the ANOVA indicates a significant treatment effect, post-hoc tests did not identify specific pairs with significant differences. "
else:
    markdown += "The analysis did not detect a significant effect of treatment on CFU counts. "

if anova_table.loc['C(Position)', 'PR(>F)'] < 0.05:
    significant_positions = posthoc_position[posthoc_position['p-tukey'] < 0.05]
    if not significant_positions.empty:
        markdown += "There is a significant effect of position on CFU counts. "
        for _, row in significant_positions.iterrows():
            markdown += f"The {row['A']} and {row['B']} positions differ significantly (p = {row['p-tukey']:.4f}). "
    else:
        markdown += "Although the ANOVA indicates a significant position effect, post-hoc tests did not identify specific pairs with significant differences. "
else:
    markdown += "The position on the slope did not significantly affect CFU counts. "

if anova_table.loc['C(Treatment):C(Position)', 'PR(>F)'] < 0.05:
    markdown += "The significant interaction between treatment and position indicates that the effect of treatments varies depending on the position on the slope."

# Calculate overall efficacy for conclusion
botector_efficacy = summary_stats[summary_stats['Treatment'] == 'Botector']['mean'].mean()
pb_efficacy = summary_stats[summary_stats['Treatment'] == 'Potassium Bicarbonate']['mean'].mean()
control_efficacy = summary_stats[summary_stats['Treatment'] == 'Control']['mean'].mean()

bot_percent = ((control_efficacy - botector_efficacy) / control_efficacy) * 100
pb_percent = ((control_efficacy - pb_efficacy) / control_efficacy) * 100

markdown += f"\n\nOverall, Botector reduced CFU counts by {bot_percent:.4f}% and Potassium Bicarbonate by {pb_percent:.4f}% compared to the Control treatment."

# Save markdown to a file
with open('cfu_analysis_report.md', 'w') as f:
    f.write(markdown)

print("Enhanced analysis complete. Check the generated plots and the markdown report.")
