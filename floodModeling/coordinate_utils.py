# coordinate_utils.py
"""
Shared utility for consistent coordinate formatting across all modules.
"""

def format_coordinate(value: float) -> str:
    """
    Format coordinate to exactly 4 decimal places with trailing zeros preserved.
    
    Parameters:
    -----------
    value : float
        Coordinate value (latitude or longitude)
        
    Returns:
    --------
    str
        Formatted string with exactly 4 decimal places
        
    Examples:
    ---------
    format_coordinate(40.8) -> "40.8000"
    format_coordinate(-73.8979) -> "-73.8979"
    format_coordinate(40.71285) -> "40.7129"
    """
    # Round to 4 decimal places
    rounded = round(float(value), 4)
    
    # Format with exactly 4 decimal places
    formatted = f"{rounded:.4f}"
    
    return formatted
