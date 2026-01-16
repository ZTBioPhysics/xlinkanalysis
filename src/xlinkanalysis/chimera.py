"""
ChimeraX integration for crosslink visualization.

This module generates commands for visualizing crosslinks in UCSF ChimeraX,
including distance measurements and pseudobond creation.
"""

import pandas as pd
from typing import Optional, List


def format_distance_commands(df: pd.DataFrame,
                             model: int = 1,
                             atom: str = 'CA') -> str:
    """
    Generate ChimeraX distance measurement commands for crosslinks.

    Args:
        df: DataFrame with 'Residue1' and 'Residue2' columns.
        model: Model number in ChimeraX session.
        atom: Atom name to measure from (typically 'CA' or 'NZ' for lysines).

    Returns:
        Multi-line string of ChimeraX distance commands.

    Example:
        >>> commands = format_distance_commands(df, model=1)
        >>> print(commands)
        distance #1:100@CA #1:150@CA
        distance #1:200@CA #1:250@CA

    Usage in ChimeraX:
        1. Copy output to clipboard
        2. In ChimeraX command line, paste commands
        3. Optionally run: distance style color green radius 0.5 dashes 0
    """
    if 'Residue1' not in df.columns or 'Residue2' not in df.columns:
        raise ValueError("DataFrame must include 'Residue1' and 'Residue2' columns.")

    commands = []
    for _, row in df.iterrows():
        res1 = int(row['Residue1'])
        res2 = int(row['Residue2'])
        cmd = f"distance #{model}:{res1}@{atom} #{model}:{res2}@{atom}"
        commands.append(cmd)

    return '\n'.join(commands)


def format_pseudobond_commands(df: pd.DataFrame,
                               model: int = 1,
                               atom: str = 'CA',
                               color: str = 'green',
                               radius: float = 0.5,
                               group_name: str = 'crosslinks') -> str:
    """
    Generate ChimeraX commands to create pseudobonds for crosslinks.

    Pseudobonds are more efficient than distance objects for large numbers
    of crosslinks and can be styled as a group.

    Args:
        df: DataFrame with 'Residue1' and 'Residue2' columns.
        model: Model number in ChimeraX session.
        atom: Atom name to connect.
        color: Pseudobond color.
        radius: Pseudobond radius.
        group_name: Name for the pseudobond group.

    Returns:
        Multi-line string of ChimeraX commands.

    Example:
        >>> commands = format_pseudobond_commands(df)
        >>> # Save to file and run in ChimeraX: open commands.cxc
    """
    commands = [f"# Crosslink pseudobonds - {len(df)} links"]
    commands.append(f"pbgroup create {group_name}")

    for _, row in df.iterrows():
        res1 = int(row['Residue1'])
        res2 = int(row['Residue2'])
        cmd = f"pbond #{model}:{res1}@{atom} #{model}:{res2}@{atom} group {group_name}"
        commands.append(cmd)

    # Style commands
    commands.append(f"setattr pbgroup:{group_name} color {color}")
    commands.append(f"setattr pbgroup:{group_name} radius {radius}")

    return '\n'.join(commands)


def format_color_by_distance(df: pd.DataFrame,
                             model: int = 1,
                             satisfied_color: str = 'green',
                             violated_color: str = 'red',
                             threshold: float = 20.0) -> str:
    """
    Generate commands to color crosslinks by satisfaction status.

    Args:
        df: DataFrame with 'Residue1', 'Residue2', and 'CA Distance' columns.
        model: Model number in ChimeraX session.
        satisfied_color: Color for crosslinks within threshold.
        violated_color: Color for crosslinks exceeding threshold.
        threshold: Distance threshold in Angstroms.

    Returns:
        Multi-line string of ChimeraX commands.
    """
    if 'CA Distance' not in df.columns:
        raise ValueError("DataFrame must include 'CA Distance' column.")

    satisfied = df[df['CA Distance'] <= threshold]
    violated = df[df['CA Distance'] > threshold]

    commands = [f"# Satisfied crosslinks ({len(satisfied)} links, <= {threshold} A)"]
    commands.append(f"pbgroup create satisfied")
    for _, row in satisfied.iterrows():
        res1, res2 = int(row['Residue1']), int(row['Residue2'])
        commands.append(f"pbond #{model}:{res1}@CA #{model}:{res2}@CA group satisfied")
    commands.append(f"setattr pbgroup:satisfied color {satisfied_color}")
    commands.append(f"setattr pbgroup:satisfied radius 0.5")

    commands.append(f"\n# Violated crosslinks ({len(violated)} links, > {threshold} A)")
    commands.append(f"pbgroup create violated")
    for _, row in violated.iterrows():
        res1, res2 = int(row['Residue1']), int(row['Residue2'])
        commands.append(f"pbond #{model}:{res1}@CA #{model}:{res2}@CA group violated")
    commands.append(f"setattr pbgroup:violated color {violated_color}")
    commands.append(f"setattr pbgroup:violated radius 0.5")

    return '\n'.join(commands)


def format_selection_command(residues: List[int], model: int = 1) -> str:
    """
    Generate a ChimeraX selection command for a list of residues.

    Args:
        residues: List of residue numbers.
        model: Model number.

    Returns:
        ChimeraX selection command string.

    Example:
        >>> cmd = format_selection_command([100, 150, 200])
        >>> print(cmd)
        select #1:100,150,200
    """
    residue_str = ','.join(str(r) for r in sorted(residues))
    return f"select #{model}:{residue_str}"


def get_crosslinked_residues(df: pd.DataFrame) -> List[int]:
    """
    Get a sorted list of all residues involved in crosslinks.

    Args:
        df: DataFrame with 'Residue1' and 'Residue2' columns.

    Returns:
        Sorted list of unique residue numbers.
    """
    residues = set(df['Residue1']).union(set(df['Residue2']))
    return sorted(residues)


def save_chimera_script(commands: str, filepath: str) -> None:
    """
    Save ChimeraX commands to a .cxc script file.

    Args:
        commands: Command string (from format_* functions).
        filepath: Output file path (should end in .cxc).

    Usage:
        In ChimeraX: open crosslinks.cxc
    """
    with open(filepath, 'w') as f:
        f.write(commands)
    print(f"ChimeraX script saved to: {filepath}")


def generate_visualization_script(df: pd.DataFrame,
                                  pdb_path: str,
                                  model: int = 1,
                                  color_by_satisfaction: bool = True,
                                  threshold: float = 20.0) -> str:
    """
    Generate a complete ChimeraX visualization script.

    Args:
        df: DataFrame with crosslink data.
        pdb_path: Path to the PDB file.
        model: Model number to assign.
        color_by_satisfaction: If True, color by distance satisfaction.
        threshold: Distance threshold for coloring.

    Returns:
        Complete ChimeraX script as a string.
    """
    commands = [
        "# Crosslink Visualization Script",
        "# Generated by xlinkanalysis",
        "",
        f"open {pdb_path}",
        "hide atoms",
        "show cartoons",
        "color bychain",
        "",
    ]

    if color_by_satisfaction and 'CA Distance' in df.columns:
        commands.append(format_color_by_distance(df, model, threshold=threshold))
    else:
        commands.append(format_pseudobond_commands(df, model))

    commands.extend([
        "",
        "# Highlight crosslinked residues",
        format_selection_command(get_crosslinked_residues(df), model),
        "color sel gold target a",
        "show sel atoms",
        "~select",
        "",
        "# Adjust view",
        "view",
    ])

    return '\n'.join(commands)
