def stringify(obj):
    """Converts all objects in a dict to str, recursively."""
    if isinstance(obj, dict):
        return {str(key): stringify(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [stringify(value) for value in obj]
    else:
        return str(obj)
