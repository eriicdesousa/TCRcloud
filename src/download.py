import json
import requests

import airr

def testserver(data):
    # Find which server has the matching rearragements
    query = {
    'filters':{
         'op':'=',
        'content': {
            'field':'repertoire_id',
            'value':data['Repertoire'][0]['repertoire_id']
                    }
            }
        }
    repositories = ['https://vdjserver.org/airr/v1',
                    'https://ipa1.ireceptor.org/airr/v1/',
                    'https://ipa2.ireceptor.org/airr/v1/',
                    'https://ipa3.ireceptor.org/airr/v1/',
                    'https://ipa4.ireceptor.org/airr/v1/',
                    'http://ipa5.ireceptor.org/airr/v1/',
                    'http://covid19-1.ireceptor.org/airr/v1/',
                    'http://covid19-2.ireceptor.org/airr/v1/',
                    'http://covid19-3.ireceptor.org/airr/v1/',
                    'http://covid19-4.ireceptor.org/airr/v1/',
                    'https://scireptor.dkfz.de/airr/v1/',
                    'http://airr-seq.vdjbase.org/airr/v1/']

    for i in repositories:
        test_url = i
        resp = requests.post(test_url + '/rearrangement', json = query)
        if str(resp) == '<Response [200]>':
            if len(resp.json()['Rearrangement']) > 0:
                host_url = i
    return host_url

def airrdownload(args):
    repertoire_file = args.repertoire
    rearrangements_file = repertoire_file[:-4]+'rearrangements.tsv'
    data = airr.load_repertoire(args.repertoire)
    repertoires = data['Repertoire']
    host_url = testserver(data)
    print(host_url)
    # Print out some Info
    print('       Info: ' + data['Info']['title'])
    print('    version: ' + str(data['Info']['version']))
    print('description: ' + data['Info']['description'])
    print('Found ' + str(len(data['Repertoire'])) + ' repertoires in \
metadata file.')

    # Query the rearrangement endpoint
    # Define a generic query object, and we will replace the repertoire_id
    # within the loop. We also only request productive rearrangements as
    # an additional filter.

    query = {
        'filters':{
            'op':'and',
            'content': [
                {
                    'op':'=',
                    'content': {
                        'field':'repertoire_id',
                        'value':'XXX'
                    }
          },
          {
                    'op':'=',
                    'content': {
                        'field':'productive',
                        'value':True
                    }
          }
      ]
        },
        'size':1000,
        'from':0
    }

    # Loop through each repertoire and query rearrangement data for
    # each. We download in chunks of 10000 because of the server 
    # limitations using the from and size parameters.

    first = True
    for r in repertoires:
        print('Retrieving rearrangements for repertoire: ' + r['repertoire_id'])
        query['filters']['content'][0]['content']['value'] = r['repertoire_id']
        query['size'] = 1000
        query['from'] = 0

        cnt = 0
        while True:
            # send the request
            resp = requests.post(host_url + '/rearrangement', json = query)
            data = resp.json()
            rearrangements = data['Rearrangement']

            # Open a file for writing the rearrangements. We do this here
            # because we need to know the full set of fields being
            # returned from the data repository, otherwise by default only
            # the required fields will be written to the file.
            if first:
                out_file = airr.create_rearrangement(rearrangements_file, fields=rearrangements[0].keys())
                first = False
            
            # save the rearrangements to a file
            for row in rearrangements:
                out_file.write(row)

            # stop when downloaded at most 10,000 rearrangements or if the
            # response doesn't return the full amount, which indicates no more
            # data. If you wanted to download all rearrangements, keep
            # looping until zero rearrangements are returned from the query.
            cnt += len(rearrangements)

            if len(rearrangements) < 1000:
               break

            # Need to update the from parameter to get the next chunk
            query['from'] = cnt

        print('Retrieved ' + str(cnt) + ' rearrangements for repertoire: ' + r['repertoire_id'])




























    # #
    # # Query the repertoire endpoint
    # #

    # # POST data is sent with the query. Here we construct an object for
    # # the query ((study_id == 'PRJNA300878') AND (locus == 'TRB'))
    
    # # Send the query
    # resp = requests.post(host_url + '/repertoire', json = query)
    
    # # The data is returned as JSON, use AIRR library to write out data
    # data = resp.json()
    # airr.write_repertoire(repertoire_file, data['Repertoire'], info=data['Info'])
    # repertoires = data['Repertoire']

    # # Print out some Info
    # print('       Info: ' + data['Info']['title'])
    # print('    version: ' + str(data['Info']['version']))
    # print('description: ' + data['Info']['description'])

    # # Save repertoires to a file
    # print('\nReceived ' + str(len(data['Repertoire'])) + ' repertoires.\n')

    # #
    # # Query the rearrangement endpoint
    # #

    # # Define a generic query object, and we will replace the repertoire_id
    # # within the loop. We also only request productive rearrangements as
    # # an additional filter.

    # if chain == 'TRA':
    #   query = {
    #       'filters':{
    #           'op':'and',
    #           'content': [
    #               {
    #                   'op':'=',
    #                   'content': {
    #                       'field':'repertoire_id',
    #                       'value':'XXX'
    #                   }
    #               },
    #               {
    #                   'op':'=',
    #                   'content': {
    #                       'field':'productive',
    #                       'value':True
    #                   }
    #               },
    #               {
    #                   'op':'=',
    #                   'content': {
    #                       'field':'locus',
    #                       'value':'TRA'
    #                   }
    #               }
    #           ]
    #       },
    #       'size':1000,
    #       'from':0
    #   }
    # elif chain == 'TRB':
    #   query = {
    #       'filters':{
    #           'op':'and',
    #           'content': [
    #               {
    #                   'op':'=',
    #                   'content': {
    #                       'field':'repertoire_id',
    #                       'value':'XXX'
    #                   }
    #               },
    #               {
    #                   'op':'=',
    #                   'content': {
    #                       'field':'productive',
    #                       'value':True
    #                   }
    #               },
    #               {
    #                   'op':'=',
    #                   'content': {
    #                       'field':'locus',
    #                       'value':'TRB'
    #                   }
    #               }
    #           ]
    #       },
    #       'size':1000,
    #       'from':0
    #   }
    # else:
    #   query = {
    #       'filters':{
    #           'op':'and',
    #           'content': [
    #               {
    #                   'op':'=',
    #                   'content': {
    #                       'field':'repertoire_id',
    #                       'value':'XXX'
    #                   }
    #               },
    #               {
    #                   'op':'=',
    #                   'content': {
    #                       'field':'productive',
    #                       'value':True
    #                   }
    #               }
    #           ]
    #       },
    #       'size':1000,
    #       'from':0
    #   }

    # # Loop through each repertoire and query rearrangement data for
    # # each. In this example, we only download 20000 rearrangements for each
    # # repertoire in chunks of 1000 using the from and size parameters.


    # first = True
    # for r in repertoires:
    #   print('Retrieving rearrangements for repertoire: ' + r['repertoire_id'])
    #   query['filters']['content'][0]['content']['value'] = r['repertoire_id']
    #   query['size'] = 1000
    #   query['from'] = 0

    #   cnt = 0
    #   while True:
    #       # send the request
    #       resp = requests.post(host_url + '/rearrangement', json = query)
    #       data = resp.json()
    #       rearrangements = data['Rearrangement']

    #       # Open a file for writing the rearrangements. We do this here
    #       # because we need to know the full set of fields being
    #       # returned from the data repository, otherwise by default only
    #       # the required fields will be written to the file.
    #       if first:
    #           out_file = airr.create_rearrangement(rearrangement_file, fields=rearrangements[0].keys())
    #           first = False
    #       # save the rearrangements to a file
    #       for row in rearrangements:
    #           out_file.write(row)

    #       # stop when downloaded at most 40,000 rearrangements or if the
    #       # response doesn't return the full amount, which indicates no more
    #       # data. If you wanted to download all rearrangements, keep
    #       # looping until zero rearrangements are returned from the query.
    #       cnt += len(rearrangements)
    #       if cnt >= 40000 or len(rearrangements) < 1000:
    #           break

    #       # Need to update the from parameter to get the next chunk
    #       query['from'] = cnt

    #   print('Retrieved ' + str(cnt) + ' rearrangements for repertoire: ' + r['repertoire_id'] + '\n')