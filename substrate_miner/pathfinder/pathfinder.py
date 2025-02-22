################################################################################# 
# PathFinder: A tool to find the pathological/molecular path for a substrate.   #
#             This tool uses the KEGG API to retrieve                           #    
#             pathways and diseases associated                                  #
#################################################################################

import os
import sys
import click
import requests

from Bio.KEGG import REST
from Bio.KEGG import Enzyme

def validate_extension(ctx, param, value):
    if value and not value.endswith(('.txt', 'text', '.csv')):
        raise click.BadParameter('File must have a .txt, .text or .csv extension')
    return value

def validate_orgs(ctx, param, value):
    if not ctx.params.get('api') and not value:
        raise click.BadParameter('Organism code is required if --api is not given.')
    return value

@click.command('pathfinder', help='Find the pathological/molecular path for a substrate.')
@click.option('--input', '-i', callback=validate_extension,\
              type=click.Path(exists=True, readable=True, resolve_path=True),
              help='Input file path')
@click.option('--output', '-o', help='Output file path')
@click.option('--api', '-a', is_flag=True, help='Use KEGG API to retrieve pathways and diseases')
@click.option('--uniprots', '-u',\
    help='UniProt ID for a protein, comma-separated for multiple IDs (e.g., P12345,Q67890) or space-separated for multiple IDs (e.g., "P12345 Q67890")')
@click.option('--orgs', '-g', default='hsa', callback=validate_orgs,\
    help='Organism code for the KEGG API (default: hsa)')
@click.pass_context
def pathfinder_main(ctx, input, output, api, uniprots, orgs):
    """
    Interface for path/biological pathway finder.
    """
    if not input and not uniprots:
        click.echo("Error: Please provide either an input file or a UniProt ID.")
        ctx.exit(1)

    if input:
        with open(input, 'r') as file:
            lines = file.readlines()
            uniprot_ids = []
            orgs = []
            for line in lines:
                parts = line.strip().split(',')
                if len(parts) == 1:
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        uniprot_ids.append(parts[0].strip())
                        orgs.append(parts[1].strip())
                    else:
                        # default organism code hsa - human
                        uniprot_ids.append(parts[0].strip())
                        orgs.append('hsa')
                elif len(parts) == 2:
                    uniprot_ids.append(parts[0].strip())
                    orgs.append(parts[1].strip())
                else:
                    click.echo("Error: Invalid input format of UniProt ID, organism code.")
                    ctx.exit(1)    
                
                # debug block
                #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                #print(uniprot_ids)
                #print(orgs)
                #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    if uniprots:
        uniprots = uniprots.replace(",", " ").split()
    
    # sanity check for output file
    if output:
        # Check if the output file exists
        if os.path.exists(output):
            os.remove(output)  # Delete the existing file
            with open(output, 'w') as file:
                file.write("UniProt_ID\tPathways\tDiseases\n")
        else:
            # Create the output file
            with open(output, 'w') as file:
                file.write("UniProt_ID\tPathways\tDiseases\n")

    if api:
        if input:
            if not output:
                output = None

            workflow_api(uniprot_ids, output)
                    
        else:
            if not output:
                output = None

            workflow_api(uniprots, output)

    else:
        ## TODO: Implement BioPython KEGG API Fix
        # Use BioPython KEGG API to retrieve pathways and diseases
        if input:
            if not output:
                output = None
                
            workflow_biokegg(uniprot_ids, orgs, output)
        else:
            if not output:
                output = None

            workflow_biokegg(uniprots, orgs, output)
            
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print(f"Pathfinder module for {uniprots} completed.")
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

def workflow_biokegg(uniprots, orgs, output) -> None:
    """
    Workflow for the BioPython KEGG API mode.
    
    Parameters:
        uniports (list): List of UniProt IDs.
        output (str): Output file path.

    Returns:
        None
    """
    
    # Retrieve pathways and diseases for each UniProt ID
    pathways = []
    diseases = []
    
    index_counter = 0
    
    for uniprot in uniprots:
        pathways = get_kegg_pathways(uniprot, orgs[index_counter])
        diseases = get_kegg_diseases(uniprot, orgs[index_counter])
        
        ## debug block
        #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        #print(pathways)
        #print(diseases)
        #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

        # stats block
        print("++++++++++++++++++++++++++++++++++++++++++++")
        print(f"UniProt ID: {uniprot}, total pathways: {len(pathways)}, total diseases: {len(diseases)}")
        print("++++++++++++++++++++++++++++++++++++++++++++")   

        # Write output to file or print to console
        if output != None:
            write_output_table(output, uniprot, pathways, diseases)
        else:
            ## Implement output to console
            click.echo(f"{uniprot}\t{'; '.join(pathways)}\t{'; '.join(diseases)}")

        index_counter += 1

def workflow_api(uniports, output) -> None:
    """
    Workflow for the API mode.
    
    Parameters:
        uniports (list): List of UniProt IDs.
        output (str): Output file path.

    Returns:
        None
    """

    uniprot_ids = uniports

    for uniprot in uniprot_ids:
        # Initialize pathways and diseases
        pathways = []
        diseases = []

        # Retrieve pathways and diseases for each UniProt ID
        pathway = get_kegg_pathways_api(uniprot)
        if "Error" in pathway:
            click.echo(f"No pathways found for UniProt ID: {uniprot}")

        disease = get_kegg_diseases_api(uniprot)
        if "Error" in disease:
            click.echo(f"No diseases found for UniProt ID: {uniprot}")

        if "Error" in pathway and "Error" in disease:
            continue

        pathways.extend(pathway)
        diseases.extend(disease)

        # stats block
        print("++++++++++++++++++++++++++++++++++++++++++++")
        print(f"UniProt ID: {uniprot}, total pathways: {len(pathways)}, "
              f"total diseases: {len(diseases)}")
        print("++++++++++++++++++++++++++++++++++++++++++++")

        # Write output to file or print to console
        if output is not None:
            write_output_table(output, uniprot, pathways, diseases)
        else:
            ## Implement output to console
            click.echo(f"{uniprot}\t{'; '.join(pathways)}\t{'; '.join(diseases)}")

def write_output_table(output, uniprot, pathways, diseases):
    """
    Write the output to a file in a tabulated format.
    Parameters:
        output (str): Output file path.
        uniprot_id (str): UniProt ID for a protein.
        pathways (list): List of KEGG pathways associated with the UniProt ID.
        diseases (list): List of KEGG diseases associated with the UniProt ID.
    """
    with open(output, 'a') as file:
        file.write(f"{uniprot}\t{'; '.join(pathways)}\t{'; '.join(diseases)}\n")

def get_kegg_pathways_api(uniprot_id) -> list:
    '''
    Retrieve KEGG pathways associated with a UniProt ID using the KEGG API.
    
    Parameters:
        uniprot_id (str): UniProt ID for a protein.

    Returns:
        list: List of KEGG pathways associated with the UniProt ID.

    Raises:
        Exception: If unable to retrieve pathways.
    '''
    # Convert UniProt ID to KEGG Gene ID
    url = f"http://rest.kegg.jp/conv/genes/uniprot:{uniprot_id}"
    response = requests.get(url)
    
    if response.status_code != 200:
        return f"Error: Unable to convert UniProt ID to KEGG gene ID. Status code {response.status_code}"
        sys.exit(1)
    
    kegg_gene_id = response.text.split("\t")[-1].strip()
    
    # Retrieve pathways associated with the KEGG gene ID
    pathway_url = f"http://rest.kegg.jp/link/pathway/{kegg_gene_id}"
    response = requests.get(pathway_url)
    
    if response.status_code != 200:
        return f"Error: Unable to retrieve pathways. Status code {response.status_code}"
        sys.exit(1)
    
    pathways = []
    if response.text.strip() != "":
        pathways = response.text.strip().split("\n")
    
    # Format pathways information
    pathway_list = []
    
    if pathways:        
        for pathway in pathways:
            pathway_id = pathway.split("\t")[1]
            pathway_list.append(pathway_id)
    
    return pathway_list

def get_kegg_diseases_api(uniprot_id) -> list:
    '''
    Retrieve KEGG diseases associated with a UniProt ID using the KEGG API.

    Parameters:
        uniprot_id (str): UniProt ID for a protein.

    Returns:
        list: List of KEGG diseases associated with the UniProt ID.
        
    Raises:
        Exception: If unable to retrieve diseases.
    '''
    # Convert UniProt ID to KEGG Gene ID
    url = f"http://rest.kegg.jp/conv/genes/uniprot:{uniprot_id}"
    response = requests.get(url)
    
    if response.status_code != 200:
        return f"Error: Unable to convert UniProt ID to KEGG gene ID. Status code {response.status_code}"
    
    kegg_gene_id = response.text.split("\t")[-1].strip()
    
    # Retrieve diseases associated with the KEGG gene ID
    disease_url = f"http://rest.kegg.jp/link/disease/{kegg_gene_id}"
    response = requests.get(disease_url)
    
    if response.status_code != 200:
        return f"Error: Unable to retrieve diseases. Status code {response.status_code}"
        sys.exit(1)
    
    diseases = []
    if response.text.strip() != "":
        diseases = response.text.strip().split("\n")
    
    # Format diseases information
    disease_list = []
    if diseases:
        for disease in diseases:
            disease_id = disease.split("\t")[1]
            disease_list.append(disease_id)
    
    # debug block
    #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    #print(disease_list)
    #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    
    return disease_list

##TODO: pending Biopython KEGG API implementation fix
def get_kegg_pathways(uniprot_id, org) -> list:
    '''
    Retrieve KEGG pathways associated with a UniProt ID using the BioPython KEGG API.
    
    Args:
        uniprot_id (str): UniProt ID for a protein.
        
    Returns:
        list: List of KEGG pathways associated with the UniProt ID.
    
    Raises:
        Exception: If unable to retrieve pathways.
    '''
    # Convert UniProt ID to KEGG Gene ID
    # TODO: implementation of orgs
    conversion_request = REST.kegg_conv(f"{org}", f"uniprot:{uniprot_id}")
    try:
        conversion_result = conversion_request.read()
    except Exception as e:
        raise(f"Error: Unable to convert UniProt ID to KEGG gene ID. {e}")
    
    if not conversion_result:
        return f"Error: Unable to convert UniProt ID to KEGG gene ID."
    
    ## debug block
    #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    #print(conversion_result)
    #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    kegg_gene_id = conversion_result.split("\t")[1].strip()
    
    ## debug block
    #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    #print(kegg_gene_id)
    #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    # Retrieve pathways associated with the KEGG gene ID
    try:
        pathway_request = REST.kegg_link("pathway", kegg_gene_id)
    except Exception as e:
        raise(f"Error: Unable to retrieve pathways. {e}")
    
    pathway_result = pathway_request.read()
    
    ## debug block
    #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    #print(pathway_result)
    #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    
    if not pathway_result:
        return f"Error: Unable to retrieve pathways."
    
    pathways = []
    if pathway_result.strip() != "":
        # Parsing the pathways
        pathways = [line.split("\t")[1] for line in pathway_result.strip().split("\n")]
        
    return pathways

def get_kegg_diseases(uniprot_id, org) -> list:
    '''
    Retrieve KEGG diseases associated with a UniProt ID using the BioPython KEGG API.

    Args:
        uniprot_id (str): UniProt ID for a protein.

    Returns:
        list: List of KEGG diseases associated with the UniProt ID.

    Raises:
        Exception: If unable to retrieve diseases.
    '''
    # Convert UniProt ID to KEGG Gene ID
    # TODO: implementation of orgs
    conversion_request = REST.kegg_conv(f"{org}", f"uniprot:{uniprot_id}")
    conversion_result = conversion_request.read()

    if not conversion_result:
        return f"Error: Unable to convert UniProt ID to KEGG gene ID."
    
    kegg_gene_id = conversion_result.split("\t")[1].strip()

    # Retrieve diseases associated with the KEGG gene ID
    try:
        disease_request = REST.kegg_link("disease", kegg_gene_id)
    except Exception as e:
        raise(f"Error: Unable to retrieve diseases. {e}")
        
    disease_result = disease_request.read()
    if not disease_result:
        return f"Error: Unable to retrieve diseases."
    
    diseases = []
    if disease_result.strip() != "":
        # Parsing the diseases
        diseases = [line.split("\t")[1] for line in disease_result.strip().split("\n")]
    
    return diseases


if __name__ == '__main__':
    pathfinder_main()