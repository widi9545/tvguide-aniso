class MapOptions:
    
    default_years = [2000, 2019]
    min_year = 2000
    max_year = 2020
    yearList = [int (x) for x in range(min_year, max_year+1)]
    columnNames = ["ID","Time","E","M","U","Loc","City","State","Lat","Lon","Issuance","Remark","Source"]
    mapType = ""
    mapTypeSelectOptions = [{'label': 'Default', 'value':'carto-positron'} ,{'label': 'Basic', 'value':'basic'}, {'label': 'Light', 'value':'light'}, {'label': 'Satellite', 'value':'satellite'}, {'label': 'Streets', 'value':'streets'}, {'label': 'Dark', 'value':'dark'}]
    tensor_addition_list_allowed = [x for x in range(97,122)]   