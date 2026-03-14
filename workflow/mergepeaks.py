import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

def translate_summit_position(row):
    starts = row[3].split(',')
    summits = row[9].split(',')
    result = [str(int(x)+int(y)) for x, y in zip(starts,summits)] 
    return ','.join(result) 

def extract_max_position(value_list):
    """Extracts the position of the maximum value from a comma-separated string."""
    values = list(map(float, value_list.split(',')))
    return values.index(max(values))

def get_corresponding_value(row):
    """Gets the value from the 8th column at the same position as the max value from the 7th column."""
    max_pos = extract_max_position(row[8])
    values = row[10].split(',')
    return values[max_pos] if max_pos < len(values) else '.'

def replace_with_max_value(value):
    """Replaces a comma-separated string with its maximum value."""
    return max(map(float, value.split(','))) if isinstance(value, str) else value

collapsed_df = pd.read_csv(args.input, sep='\t', header=None)

# Replace the 6th column with its maximum value and get the corresponding value
collapsed_df[10] = collapsed_df.apply(translate_summit_position, axis=1)
collapsed_df[11] = collapsed_df.apply(get_corresponding_value, axis=1)
collapsed_df[7] = collapsed_df[7].apply(replace_with_max_value)
collapsed_df[8] = collapsed_df[8].apply(replace_with_max_value)

collapsed_df.drop(collapsed_df.columns[3], axis=1, inplace=True)
collapsed_df.drop(collapsed_df.columns[3], axis=1, inplace=True)

# Insert a new 4th column with unique macs2_peak values
collapsed_df.insert(3, 'name', [f"macs2_peak_{i+1}" for i in range(len(collapsed_df))])

# Insert a new column at the 6th position with '.'
collapsed_df.insert(5, '.', '.')

# calculate new summit for 9th column
collapsed_df[9] = collapsed_df[11].astype('int')-collapsed_df[1]
collapsed_df.drop(collapsed_df.columns[10], axis=1, inplace=True)
collapsed_df.drop(collapsed_df.columns[10], axis=1, inplace=True)


# Write output to file
collapsed_df.to_csv(args.output, sep='\t', index=False, header=False)



