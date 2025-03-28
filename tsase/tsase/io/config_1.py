import configparser
import yaml
import os, sys
import copy
from icecream import ic

class ConfigClass:
    """
    Parameters are stored in a dictionary:
      {'run_type': 'trainFF', 
       'trajectory_file': 'train.traj',
       'energy_training': True,
       'force_training': True,
       'fp_type': 'BP',
       'parameter_file': 'fpParas.dat',
       'model_type': 'neural_network',
       'loss_type': 'rmse',
       'energy_coefficient': 1.0,
       'force_coefficient': 0.02,
       'optimizer_type': 'LBFGS',
       'loss_convergence': 0.01,
       'fp_paras': <__main__.FingerprintsParas object at 0x7f1c79867fd0>}

        To fetch the particular parameter, use 'config.config[<nameOFparameter>]'
           i.e., fetch run_type:
            1. Define ConfigClass object and read 'config.yaml' (default values)
               config = ConfigClass()
            2. read 'config.ini' (user defined) 
               config.initialize()
            3. Fetch the parameter
               config.config['run_type']
        Fingerprint parameters are sorted as an object FingerprintsParas. It can be fetched with
           1. Fetch the FingerprintsParas object
              fp = config.config['fp_paras']
           2. Fetch the dictionary containing the fingerprint object
              fp.fp_paras
              it will give:
                {'H': [G1, G1, ..., G2, ...],
                 'Pd': [G1, G1, ..., G2, ...]
                }
              where G1 and G2 are objects containing the necessary parameters.
               i.e., neighbor, eta, Rs for G1
                     neighbor1, neighbor2, eta, zeta, lambda, thetas for G2

    """
    def __init__(self):
        self.init_done = False
        self.cwd = os.getcwd() + '/'
        self.config = {}

        self.parser = configparser.ConfigParser()

        yaml_file = open(os.path.join(os.path.dirname(__file__), 'config.yaml'))
        self.config_defaults = yaml.load(yaml_file, Loader=yaml.BaseLoader)
        ic("self.config_defaults: ",self.config_defaults)
        yaml_file.close()

        # Pass default settings to parser
        for key in self.config_defaults:
            kattr = self.config_defaults[key]
            self.parser.set('DEFAULT', key.lower(), kattr['default'])
                #if 'values' in kattr:
                #    for value in kattr['values']:
                #        ck.values.append(value)
            #self.section_options[sectionName] = copy.deepcopy(self.parse.options(sectionName))
        #print(self.config_defaults)

    # Check legality and set values
    def setValues(self, data_type, vname, section='DEFAULT', values=None):
        config_error = False
        if data_type == 'integerlist':
            try:
                self.config[vname] = self.parser.get(section, vname)
                self.config[vname] = tuple([int(field) for field in self.config[vname].split()])
            except:
                config_error = True
                sys.stderr.write('Option "%s" in section "%s" should be an integer\n' % (vname,section))
                pass
        if data_type == 'floatlist':
            try:
                self.config[vname] = self.parser.get(section, vname)
                self.config[vname] = tuple([float(field) for field in self.config[vname].split()])
            except:
                config_error = True
                sys.stderr.write('Option "%s" in section "%s" should be a float\n' % (vname,section))
                pass
        elif data_type == 'integer':
            try:
                self.config[vname] = self.parser.getint(section, vname)
            except:
                config_error = True
                sys.stderr.write('Option "%s" in section "%s" should be a integer\n' % (vname,section))
        elif data_type == 'float':
            try:
                self.config[vname] = self.parser.getfloat(section, vname)
            except:
                config_error = True
                sys.stderr.write('Option "%s" in section "%s" should be a float\n' % (vname,section))
        elif data_type == 'boolean':
            try:
                self.config[vname] = self.parser.getboolean(section, vname)
            except:
                config_error = True
                sys.stderr.write('Option "%s" in section "%s" should be a boolean\n' % (vname,section))
                pass
        elif data_type == 'string':
            if values is not None:
                pvalue = self.parser.get(section, vname)
                if pvalue not in values:
                    config_error = True
                    sys.stderr.write('Option "%s" should be one of: %s\n' % (vname, ", ".join(values)))
                else:
                    self.config[vname] = pvalue
            else:
                self.config[vname] = self.parser.get(section, vname)
        return config_error


    def initialize(self, config_file="config.ini"):
        if os.path.isfile(config_file):
            self.parser.read(config_file)
            self.config_path = os.path.abspath(config_file)
        else:
            # TODO: config_path is not an attribute yet
            # print("Specified configuration file %s does not exist" % ''.join([self.cwd,config_file]), sys.stderr)
            print("Specified configuration file %s does not exist" % ''.join(os.path.abspath(config_file)), sys.stderr)
            sys.exit(2)

        # Check that all options in config.ini are in the configparser
        ic(self)
        ic(self.parser)
        ic(self.config.sections())
        b = self.parser.options("DEFAULT")
        section_diff = list(set(b) - set(self.config_defaults['options'].keys()))
        if len(section_diff) > 0:
            config_error = True
            sys.stderr.write('unknown option "%s" in section "%s"\n' % (", ".join(section_diff), section))

        # Check type of input settings and assignments
        defaults = self.config_defaults['options']
        for op in self.parser.options():
            if 'values' in defaults[op]:
                config_err = self.setValues(defaults[op]['kind'], op, values=defaults[op]['values'])
            else:
                config_err = self.setValues(defaults[op]['kind'], op)

        if config_err:
            sys.exit(2)

if __name__ == "__main__":
    config = ConfigClass()
    config.initialize()
    print("Parameters are stored in a dictionary:")
    print("  ",config.config)
    print("    To fetch the particular parameter, use 'config.config[<nameOFparameter>]'")
