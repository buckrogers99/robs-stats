import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
import pingouin as pg  # For post-hoc analysis

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

# Two-way ANOVA
model = ols('CFU ~ C(Treatment) + C(Position) + C(Treatment):C(Position)', data=tidy_df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("Two-Way ANOVA Results:")
print(anova_table)

# Post-hoc tests
# For Treatment
posthoc_treatment = pg.pairwise_tukey(data=tidy_df, dv='CFU', between='Treatment')
print("\nTukey's HSD Post-hoc Test for Treatment:")
print(posthoc_treatment)

# For Position
posthoc_position = pg.pairwise_tukey(data=tidy_df, dv='CFU', between='Position')
print("\nTukey's HSD Post-hoc Test for Position:")
print(posthoc_position)

# For interaction
print("\nPost-hoc for interaction (Treatment × Position):")
for pos in tidy_df['Position'].unique():
    subset = tidy_df[tidy_df['Position'] == pos]
    if len(subset) > 0:
        print(f"\nPosition: {pos}")
        posthoc = pg.pairwise_tukey(data=subset, dv='CFU', between='Treatment')
        print(posthoc)

# Create visualizations
# 1. Box plot showing treatments by position
plt.figure(figsize=(12, 6))
sns.boxplot(x='Position', y='CFU', hue='Treatment', data=tidy_df)
plt.title('CFU Counts by Position and Treatment')
plt.ylabel('Colony Forming Units (CFU)')
plt.savefig('cfu_boxplot.png')
plt.close()

# 2. Interaction plot
plt.figure(figsize=(10, 6))
interaction = tidy_df.groupby(['Position', 'Treatment'])['CFU'].mean().reset_index()
interaction_pivot = interaction.pivot(index='Position', columns='Treatment', values='CFU')
interaction_pivot.plot(marker='o', ax=plt.gca())
plt.title('Interaction Plot: Treatment × Position')
plt.ylabel('Mean CFU Count')
plt.grid(True, linestyle='--', alpha=0.7)
plt.savefig('interaction_plot.png')
plt.close()

# 3. Bar plot with error bars
plt.figure(figsize=(12, 6))
summary = tidy_df.groupby(['Treatment', 'Position'])['CFU'].agg(['mean', 'sem']).reset_index()
summary_pivot = summary.pivot(index='Treatment', columns='Position')
summary_pivot.columns = [f'{col[1]}_{col[0]}' for col in summary_pivot.columns]

bar_width = 0.25
positions = np.arange(len(summary_pivot.index))

fig, ax = plt.subplots(figsize=(12, 6))

# Plot each position as a group of bars
for i, pos in enumerate(['Top', 'Middle', 'Bottom']):
    means = summary_pivot[f'{pos}_mean']
    errors = summary_pivot[f'{pos}_sem']
    ax.bar(positions + i*bar_width - bar_width, means, width=bar_width, 
           yerr=errors, capsize=5, label=pos)

ax.set_xticks(positions)
ax.set_xticklabels(summary_pivot.index)
ax.set_ylabel('Mean CFU Count')
ax.set_title('Mean CFU Count by Treatment and Position')
ax.legend()
plt.tight_layout()
plt.savefig('cfu_barplot.png')

print("Analysis complete. Check the generated plots.")