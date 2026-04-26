#!/usr/bin/env python3
"""
Sequence Alignment Visualization Tool - Consensus-based Coloring
"""

import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import Counter

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def read_fasta(filename):
    """Read FASTA format file"""
    sequences = {}
    current_id = None
    current_seq = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def get_variable_positions(sequences):
    """Get positions with variations"""
    seq_list = list(sequences.values())
    alignment_length = len(seq_list[0])
    variable_positions = []

    for i in range(alignment_length):
        column = [seq[i] for seq in seq_list]
        if len(set(column)) > 1:
            variable_positions.append(i)

    return variable_positions


def group_consecutive_positions(positions):
    """
    Group consecutive positions into blocks
    Returns: list of (start_pos, end_pos) tuples
    """
    if not positions:
        return []

    blocks = []
    start = positions[0]
    prev = positions[0]

    for pos in positions[1:]:
        if pos != prev + 1:
            blocks.append((start, prev))
            start = pos
        prev = pos

    blocks.append((start, prev))
    return blocks


def get_column_consensus(sequences, position):
    """
    Get consensus information for a column
    Returns: (primary, secondary, primary_count, secondary_count)
    """
    seq_list = list(sequences.values())
    column = [seq[position].upper() for seq in seq_list]

    counter = Counter(column)
    most_common = counter.most_common(2)

    if len(most_common) == 0:
        return None, None, 0, 0

    primary = most_common[0][0]
    primary_count = most_common[0][1]

    if len(most_common) > 1:
        secondary = most_common[1][0]
        secondary_count = most_common[1][1]
    else:
        secondary = None
        secondary_count = 0

    return primary, secondary, primary_count, secondary_count


def get_residue_color(residue, primary, secondary):
    """
    Return color based on whether residue is primary or secondary
    Primary: Blue #4169E1 (alpha=1.0)
    Secondary: Blue #4169E1 (alpha=0.7, 70% transparency)
    Other: White #FFFFFF
    Returns: (color, alpha)
    """
    residue = residue.upper()

    if residue == primary:
        return '#4169E1', 1.0  # Blue, fully opaque
    elif residue == secondary and secondary is not None:
        return '#4169E1', 0.5  # Blue, 70% opacity
    else:
        return '#FFFFFF', 1.0  # White


def visualize_alignment(sequences, positions_to_show, residues_per_line):
    """Visualize sequence alignment with block grouping"""
    num_sequences = len(sequences)

    # Group consecutive positions into blocks
    blocks = group_consecutive_positions(positions_to_show)

    # Calculate layout
    block_gap = 3  # Gap between non-consecutive blocks (in cell units)

    line_blocks = []  # Each element: list of (block_index, x_offset) for that line
    current_line = []
    current_line_length = 0

    for block_idx, (start, end) in enumerate(blocks):
        block_length = end - start + 1

        # Check if block fits in current line
        if current_line_length > 0 and current_line_length + block_gap + block_length > residues_per_line:
            # Start new line
            line_blocks.append(current_line)
            current_line = [(block_idx, 0)]
            current_line_length = block_length
        else:
            # Add to current line
            x_offset = current_line_length + (block_gap if current_line else 0)
            current_line.append((block_idx, x_offset))
            current_line_length = x_offset + block_length

    if current_line:
        line_blocks.append(current_line)

    num_lines = len(line_blocks)

    # Calculate maximum line length
    max_line_length = max(
        [blocks_in_line[-1][1] + (blocks[blocks_in_line[-1][0]][1] - blocks[blocks_in_line[-1][0]][0] + 1)
         for blocks_in_line in line_blocks])

    # Define SQUARE cell size
    cell_size = 1.0  # One unit = one square cell

    # Calculate margins and spacing (in cell units)
    left_margin = 12  # Space for sequence names
    right_margin = 2
    top_margin = 5  # Increased space for rotated position numbers
    bottom_margin = 1
    vertical_gap = 2  # Gap between lines
    position_number_height = 3  # Increased height for rotated numbers

    # Calculate total width and height in cell units
    total_width = left_margin + max_line_length * cell_size + right_margin
    total_height = (top_margin +
                    num_lines * (num_sequences * cell_size + position_number_height + vertical_gap) +
                    bottom_margin)

    # Create figure with SQUARE aspect ratio
    # Use dpi to control resolution, figure size determines the aspect ratio
    dpi = 100
    fig_width = total_width * 0.25  # Scale factor for figure size
    fig_height = total_height * 0.25

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)

    # Set axis limits to ensure square cells
    ax.set_xlim(0, total_width)
    ax.set_ylim(0, total_height)
    ax.set_aspect('equal', adjustable='box')  # KEY: ensure equal aspect ratio
    ax.axis('off')

    # Calculate font size based on cell size
    font_size = 8  # Fixed font size for residues
    name_font_size = 7
    number_font_size = 7

    seq_ids = list(sequences.keys())
    seq_list = list(sequences.values())

    # Starting y position (from bottom)
    y_start = bottom_margin

    # Draw each line from bottom to top
    for line_idx in range(num_lines - 1, -1, -1):
        blocks_in_line = line_blocks[line_idx]

        # Calculate y position for this line
        line_offset = (num_lines - 1 - line_idx)
        y_base = y_start + line_offset * (num_sequences * cell_size + position_number_height + vertical_gap)

        # Draw blocks in this line
        for block_idx, x_offset in blocks_in_line:
            start_pos, end_pos = blocks[block_idx]
            block_positions = list(range(start_pos, end_pos + 1))

            # Draw position numbers at top of block with 90-degree rotation
            for local_idx, pos in enumerate(block_positions):
                x_pos = left_margin + (x_offset + local_idx) * cell_size

                # Show first, last, and every 5th position
                if local_idx == 0 or local_idx == len(block_positions) - 1 or (pos + 1) % 5 == 0:
                    ax.text(x_pos + cell_size / 2,
                            y_base + num_sequences * cell_size + 0.3,
                            str(pos + 1),
                            ha='center', va='bottom',
                            fontsize=number_font_size,
                            fontweight='normal',
                            color='#333333',
                            fontfamily='Arial',
                            rotation=90)  # Added 90-degree rotation

            # Draw sequences
            for seq_idx, (seq_id, seq) in enumerate(zip(seq_ids, seq_list)):
                y_pos = y_base + (num_sequences - seq_idx - 1) * cell_size

                # Draw sequence name (only for first block in line)
                if block_idx == blocks_in_line[0][0]:
                    seq_name = seq_id.split('/')[0] if '/' in seq_id else seq_id
                    # Limit name length
                    if len(seq_name) > 15:
                        seq_name = seq_name[:12] + '...'

                    ax.text(0.5, y_pos + cell_size / 2,
                            seq_name,
                            ha='left', va='center',
                            fontsize=name_font_size,
                            fontfamily='Arial',
                            color='black')

                # Draw residues
                for local_idx, pos in enumerate(block_positions):
                    residue = seq[pos].upper()
                    x_pos = left_margin + (x_offset + local_idx) * cell_size

                    # Get consensus info
                    primary, secondary, _, _ = get_column_consensus(sequences, pos)
                    bg_color, alpha = get_residue_color(residue, primary, secondary)

                    # Draw SQUARE background - NO edgecolor to remove white lines
                    rect = patches.Rectangle((x_pos, y_pos),
                                             cell_size, cell_size,
                                             linewidth=0,  # Changed from 0.5 to 0
                                             edgecolor='none',  # Changed from '#CCCCCC' to 'none'
                                             facecolor=bg_color,
                                             alpha=alpha)
                    ax.add_patch(rect)

                    # Draw residue character (always black, centered)
                    ax.text(x_pos + cell_size / 2,
                            y_pos + cell_size / 2,
                            residue,
                            ha='center', va='center',
                            fontsize=font_size,
                            fontfamily='Arial',
                            fontweight='bold',
                            color='black')

    # Add legend
    from matplotlib.patches import Patch
    legend_y = total_height - 1.5
    legend_elements = [
        Patch(facecolor='#4169E1', edgecolor='none', alpha=1.0, label='Primary Consensus'),
        Patch(facecolor='#4169E1', edgecolor='none', alpha=0.7, label='Secondary Consensus (70% opacity)'),
        Patch(facecolor='#FFFFFF', edgecolor='#CCCCCC', label='Other Residues'),
    ]
    legend = ax.legend(handles=legend_elements,
                       loc='upper right',
                       fontsize=name_font_size + 1,
                       frameon=True,
                       fancybox=False,
                       shadow=False,
                       bbox_to_anchor=(1.0, 1.0),
                       bbox_transform=ax.transData)

    # Add title
    title = f'Sequence Alignment Visualization ({num_sequences} sequences, {len(positions_to_show)} positions)'
    ax.text(total_width / 2, total_height - 1,
            title,
            ha='center', va='top',
            fontsize=name_font_size + 3,
            fontweight='bold')

    plt.tight_layout(pad=0.5)
    return fig


def main():
    parser = argparse.ArgumentParser(
        description='Sequence alignment visualization tool with consensus-based coloring',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Show all positions, 50 residues per line
  %(prog)s -i alignment.fasta -n 50

  # Show only variable positions, 30 residues per line
  %(prog)s -i alignment.fasta -diff -n 30

Color scheme:
  Blue (100%% opacity)  - Primary consensus residue (most frequent in column)
  Blue (70%% opacity)   - Secondary consensus residue (second most frequent)
  White                - Other residues

Note: 
  - Consecutive positions are grouped together without gaps
  - Non-consecutive blocks are separated by gaps
  - All text is displayed in black
  - Each residue is shown in a SQUARE cell
  - Use the matplotlib toolbar to save the figure interactively
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTA alignment file')
    parser.add_argument('-diff', '--diff-only', action='store_true',
                        help='Show only columns with variations')
    parser.add_argument('-n', '--residues-per-line', type=int, default=50,
                        help='Maximum residues per line (default: 50)')

    args = parser.parse_args()

    # Read alignment file
    print(f"Reading alignment file: {args.input}")
    try:
        sequences = read_fasta(args.input)
    except Exception as e:
        print(f"✗ Error reading file: {e}")
        return

    if not sequences:
        print("✗ No sequences found in file")
        return

    # Check sequence lengths
    seq_lengths = [len(seq) for seq in sequences.values()]
    if len(set(seq_lengths)) > 1:
        print(f"✗ Warning: Inconsistent sequence lengths! Range: {min(seq_lengths)}-{max(seq_lengths)}")
        return

    print(f"✓ Successfully read {len(sequences)} sequences of length {seq_lengths[0]}")

    # Determine positions to show
    if args.diff_only:
        positions_to_show = get_variable_positions(sequences)
        print(f"✓ Found {len(positions_to_show)} variable positions")
    else:
        positions_to_show = list(range(seq_lengths[0]))
        print(f"✓ Showing all {len(positions_to_show)} positions")

    if len(positions_to_show) == 0:
        print("✗ No positions to display!")
        return

    # Visualize
    print("Generating visualization...")
    print("Use the matplotlib toolbar to zoom, pan, and save the figure.")
    fig = visualize_alignment(sequences, positions_to_show, args.residues_per_line)
    plt.show()


if __name__ == '__main__':
    main()