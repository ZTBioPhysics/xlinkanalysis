"""
Configuration loading utilities.

This module handles loading YAML configuration files for domain definitions
and analysis parameters.
"""

from pathlib import Path
from typing import Dict, Any, Optional
import yaml


def load_yaml(filepath: str) -> Dict[str, Any]:
    """
    Load a YAML configuration file.

    Args:
        filepath: Path to the YAML file.

    Returns:
        Dictionary containing the configuration.

    Raises:
        FileNotFoundError: If the config file doesn't exist.
        yaml.YAMLError: If the file contains invalid YAML.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Configuration file not found: {filepath}")

    with open(filepath, 'r') as f:
        return yaml.safe_load(f)


def load_domains(filepath: str) -> Dict[str, Dict]:
    """
    Load domain definitions from a YAML file.

    Args:
        filepath: Path to the domains YAML file.

    Returns:
        Dictionary mapping domain names to their properties (ranges, color, etc.)

    Example:
        >>> domains = load_domains("config/domains.yaml")
        >>> domains['NTD']['ranges']
        [[1, 994]]
    """
    config = load_yaml(filepath)
    return config.get('domains', {})


def load_analysis_config(filepath: str) -> Dict[str, Any]:
    """
    Load analysis configuration from a YAML file.

    Args:
        filepath: Path to the analysis config YAML file.

    Returns:
        Dictionary containing analysis parameters.
    """
    return load_yaml(filepath)


def get_default_config_path(config_name: str) -> Path:
    """
    Get the path to a default config file in the package's config directory.

    Args:
        config_name: Name of the config file (e.g., 'domains.yaml')

    Returns:
        Path to the config file.
    """
    package_dir = Path(__file__).parent.parent.parent.parent
    return package_dir / "config" / config_name


def get_domain_ranges(domains_config: Dict[str, Dict]) -> Dict[str, list]:
    """
    Extract just the ranges from a domains configuration.

    Args:
        domains_config: Full domains configuration dictionary.

    Returns:
        Dictionary mapping domain names to their residue ranges.

    Example:
        >>> domains = load_domains("config/domains.yaml")
        >>> ranges = get_domain_ranges(domains)
        >>> ranges['NTD']
        [[1, 994]]
    """
    return {name: info['ranges'] for name, info in domains_config.items()}
