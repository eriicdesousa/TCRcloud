def jsonfile(value):
    if not value.endswith('.yml'):
        if not value.endswith('.yaml'):
            if not value.endswith('.json'):
                raise argparse.ArgumentTypeError(
                    'argument filename must be ".yaml", ".yml" or ".json"')
    return value

# print("""
#                   _____ ____ ____      _                 _ 
#                  |_   _/ ___|  _ \ ___| | ___  _   _  __| |
#                    | || |   | |_) / __| |/ _ \| | | |/ _` |
#                    | || |___|  _ < (__| | (_) | |_| | (_| |
#                    |_| \____|_| \_\___|_|\___/ \__,_|\__,_|
#    """) 