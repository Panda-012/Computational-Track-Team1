import os

# List of IDs and their corresponding cancer subtypes
id_to_subtype = {
    "GSM1116169": "Normal",
    "GSM1116170": "Normal",
    "GSM1116171": "Normal",
    "GSM1116172": "Normal",
    "GSM1116173": "Normal",
    "GSM1116174": "Normal",
    "GSM1116175": "Normal",
    "GSM1116176": "Normal",
    "GSM1116177": "Normal",
    "GSM1116178": "Normal",
    "GSM1116179": "Normal",
    "GSM1116180": "Luminal_A",
    "GSM1116181": "Luminal_A",
    "GSM1116182": "Luminal_B",
    "GSM1116183": "Luminal_B",
    "GSM1116184": "Luminal_A",
    "GSM1116185": "Luminal_B",
    "GSM1116186": "Luminal_B",
    "GSM1116187": "Luminal_A",
    "GSM1116188": "Luminal_B",
    "GSM1116189": "Luminal_B",
    "GSM1116190": "Luminal_A",
    "GSM1116191": "Luminal_A",
    "GSM1116192": "Luminal_B",
    "GSM1116193": "Luminal_B",
    "GSM1116194": "Luminal_A",
    "GSM1116195": "Luminal_A",
    "GSM1116196": "Luminal_B",
    "GSM1116197": "Luminal_B",
    "GSM1116198": "Luminal_B",
    "GSM1116199": "Luminal_A",
    "GSM1116200": "Luminal_A",
    "GSM1116201": "Luminal_B",
    "GSM1116202": "Luminal_B",
    "GSM1116203": "Luminal_A",
    "GSM1116204": "Luminal_A",
    "GSM1116205": "Luminal_A",
    "GSM1116206": "Luminal_B",
    "GSM1116207": "Luminal_B",
    "GSM1116208": "Luminal_A",
    "GSM1116209": "Luminal_A",
    "GSM1116210": "Luminal_B",
    "GSM1116211": "Luminal_B",
    "GSM1116212": "Luminal_B",
    "GSM1116213": "Luminal_A",
    "GSM1116214": "Luminal_B",
    "GSM1116215": "Luminal_B",
    "GSM1116216": "Luminal_A",
    "GSM1116217": "Luminal_A",
    "GSM1116218": "Luminal_B",
    "GSM1116219": "Luminal_B",
    "GSM1116220": "Luminal_A",
    "GSM1116221": "Luminal_A",
    "GSM1116222": "Luminal_A",
    "GSM1116223": "Luminal_A",
    "GSM1116224": "Luminal_A",
    "GSM1116225": "Luminal_A",
    "GSM1116226": "Luminal_B",
    "GSM1116227": "Luminal_B",
    "GSM1116228": "Luminal_A",
    "GSM1116229": "Luminal_B",
    "GSM1116230": "Luminal_B",
    "GSM1116231": "Luminal_A",
    "GSM1116232": "Luminal_A",
    "GSM1116233": "Luminal_B",
    "GSM1116234": "Luminal_A",
    "GSM1116235": "Luminal_A",
    "GSM1116236": "Luminal_B",
    "GSM1116237": "Luminal_B",
    "GSM1116238": "Luminal_B"
}

# Initialize a dictionary to store counts for each subtype
subtype_counts = {
    "Normal": 0,
    "Luminal_A": 0,
    "Luminal_B": 0
}

# Get all .CEL files in the directory
files = [f for f in os.listdir() if f.endswith(".CEL")]

# Rename files based on IDs and subtype
for file in files:
    file_parts = file.split('_')  # Split filename by underscores
    for id_key, subtype in id_to_subtype.items():
        if id_key == file_parts[0]:  # Check the first part of the filename
            subtype_counts[subtype] += 1
            new_name = f"{subtype}_{subtype_counts[subtype]}.CEL"  # Reconstruct filename
            os.rename(file, new_name)
            break

# Print the subtype counts
print(subtype_counts)
