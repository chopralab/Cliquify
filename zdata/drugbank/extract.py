import requests

api_key = 'fong14@purdue.edu:Minatozaki123'
api_url = 'https://api.drugbank.ca/v1/drugs'

headers = {
    'Authorization': f'Bearer {api_key}'
}

# response = requests.get(api_url, headers=headers)
response = requests.get(api_url)

import json


try:
    data = json.loads(response.text)
    # Rest of the code to extract SMILES
except json.JSONDecodeError as e:
    print("Error decoding JSON response:")

    with open('drugbank_smiles.txt', 'w') as f:
        f.write(response.text)
    print(response.text)
    
    raise e


smiles_list = [entry['structures']['canonical_smiles'] for entry in data['data']]

