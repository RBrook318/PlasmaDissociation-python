try:
    import pyqchem.tools
    print("tools imported successfully")
except ImportError as e:
    print(f"ImportError: {e}")