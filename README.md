# RNAchat
Understanding the global interactions among cellular pathways and types is essential for unraveling the complex biological mechanisms that underlie various diseases. While many existing studies have primarily focused on omics data, they often fail to integrate clinical heterogeneity, which is crucial for a comprehensive understanding of disease biology. Although significant progress has been made in researching crosstalks between inter-pathways and inter-cells, and tools like CellChat have emerged to analyze these interactions, they have notable limitations. These include the assumption that mRNA expression directly reflects protein activity, a heavy reliance on pre-established ligand-receptor interactions, and an inability to connect “chats” with clinical phenotypes. Furthermore, existing methods have limitations in distinguishing between paracrine signaling (interactions between different cell types) and autocrine signaling (self-signaling within the same cell type). These limitations hinder the ability to translate molecular insights into phenotypic outcomes, thereby limiting clinical applicability. To provide another solution for these challenges, we introduce RNAchat, a novel approach that utilizes machine learning techniques to identify metapathways for crosstalks among pathways at not only the bulk level but also cell types at the single-cell level, incorporating both clinical and multi-omics data.

RNAchat provides another solution by providing an interactive, reproducible platform designed to analyse metapathways between inter-pathways and inter-cell types using multi-omics level data in a clinical context. This tool enables researchers to integrate diverse datasets, conduct exploratory analyses, and identify cooperated pathways. The platform facilitates hypothesis generation and produces intuitive visualizations to support further investigation.
As a proof of concept, we applied RNAchat to connect omics data to pain, drug resistance in rheumatoid arthritis (RA), disease severity to RA and Heart Failure (HF). We discovered cooperated molecular pathways associated with different treatments using SHAP (Shapley Additive Explanations) values. Additionally, we demonstrate how to use it to find cellular interactions at single cell level.


## System Implementation and Package dependence
Operating system(s): Platform independent
Compatible bBrowsers: Firefox/Chrome
Programming language: Python, JavaScript
Other requirements: Python >= 3.11, Django >= 4.2, Nginx(optional)
Any restrictions to use by non-academics: licence needed
### 1. Complete code for migration:
To copy RNAcompare from GitHub and make it runnable on your local machine, you can follow these steps:
#### Step 1: Clone the Repository
First, clone the repository from GitHub to your local machine.
```bash
https://github.com/tangmingcan/RNAchat.git
cd RNAchat
```
#### Step 2: Set Up a Virtual Environment(for development, we use python 3.11)
It's a good practice to use a virtual environment to manage your project's dependencies.
```python
# Install virtualenv if you haven't already
pip install virtualenv

# Create a virtual environment
python3 -m virtualenv venv

# Activate the virtual environment
# On Windows
venv\Scripts\activate
# On macOS/Linux
source venv/bin/activate
```
#### Step 3: Install Dependencies
```bash
pip install -r requirements.txt
```
#### Step 4: Configure the Django Settings
Ensure the Django settings file is configured correctly. The default settings for SQLite should be fine if you're running it locally.

Open the settings.py file in your Django project directory and check the database settings and modify the uploaded folder:
```python
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / "db.sqlite3",
    }
}

# change to your own folder for storing user specific uploaded expression and clinic files and storing the shared cohorts data
MEDIA_URL = "/media/"
MEDIA_ROOT = os.path.join(BASE_DIR, "uploaded")

# change to your local server if it is not 127.0.0.1
ALLOWED_HOSTS = ["127.0.0.1"]
```
You might also want to change the directory for user uploaded files from the default value to some other directory that you have write access to:
```python
MEDIA_ROOT = os.path.join(BASE_DIR, "uploaded")
```
#### Step 5: Run Migrations to configure database
In order to create the database db.sqlite3 and configure the schema, run the following:
```bash
python manage.py migrate
```
#### Step 6: Create a Superuser (if needed)
If you need to create a superuser for accessing the Django admin interface, you can do so with:
```bash
python manage.py createsuperuser
```
#### Step 7: Configure Celery and Redis
The project is built based on Django + Celery + Redis to avoid some issues of fig.save for matplotlib. In order to use Celery and Redis, please have a reference about how to set up it: https://realpython.com/asynchronous-tasks-with-django-and-celery/

Install Redis server:
```bash
# On Linux/WSL2(Windows)
sudo apt install redis
# On macOS
brew install redis
```
To start the Redis server, open another console:
```bash
redis-server --port 8001
```
Then open another console for test:
```bash
redis-cli -p 8001
```
If successfully connected, you should see a cursor with 127.0.0.1:8001>.

Note: Celery settings are defined in djangoproject/settings.py and djangoproject/Celery.py, with the corresponding port number and serialization method for redis.

Open another console to start Celery:
```bash
python3 -m celery -A djangoproject worker -l info
```
#### Step 8: Run the Django Development Server
Finally, run the Django development server to verify that everything is set up correctly.
```bash
python manage.py runserver
```
Open your web browser and go to http://127.0.0.1:8000/ to see your Django project running.

##### 8.1 Role management for users
Similar to RNAcare, you can do role management for different users to control their access to shared datasets.

#### Step 9(optional): Nginx setting
Edit your Nginx configuration file (e.g., /etc/nginx/sites-available/your_site.conf) and ensure it includes the following directive inside the appropriate location or server block:
```bash
proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
proxy_set_header X-Real-IP $remote_addr;
```
Reload Nginx to apply the changes:`sudo systemctl reload nginx`
#### Step 10: Constants settings
In constants.py, there are a few global parameters:
```
ALLOW_UPLOAD = True # This is to control the platform whether to allow user to upload data
SHAP_PLOT_DATASET='X' # This is for SHAP importance plot and dependence plot based on X/X_train/X_test dataset.
```

## System Introduction
### Platform Features:
##### Metapathway Analysis
##### Associatioin Analysis
##### Double Machine Learning (DML) Analysis
##### Data Visualisation
### Materials and Methods
#### Workflow overview
##### A - Data Upload 
The user has the option to import preprocessed data from RNAcompare, or upload their own iCA pre-processed data on the fly combined with clinical data which may include cell type information for scRNA analysis. Clinical data are uploaded as a table.
##### B - Metapathway Analysis
Based on the first-level independent component analysis (ICA24) processed data, we perform a second-level ICA on the metagenes. We hypothesize that the pathways represented by metagenes interact with one another, forming higher-order structures referred to as metapathways. Following the second-level ICA, we can extract the components of these metapathways, each composed of multiple metagenes. In the case of single-cell RNA sequencing (scRNA-seq) data, individual cell types are treated as specialized metagenes during the decomposition process. Consequently, the resulting metapathways may incorporate cell type information, reflecting ICA-derived correlations between cell types and pathways within the same metapathway. If a metapathway contains only a single cell type, it is indicative of autocrine signaling; if multiple cell types are present, it suggests paracrine signaling. Conversely, metapathways devoid of cell type representation correspond to generalized pathway interactions that are not cell type-specific.
Metapathways can be directly downloaded after processing.

##### C – Data Visualisation I
In this section, users will perceive the relationships of metapathways. We provide 2 ways: Pearson Correlation Matrix25 and matrixplot scaled by metapathways/batches. When uses want to compare the expression level of different metapathways within the same tissue and omics level, the option of ‘scaled by metapathways’ is recommended; otherwise, ‘scaled by batches’ is recommended. 

##### D – Associatioin Analysis
Here, the user has the option to associate selected clinical data parameters with metapathway data in order to study phenotypes based on part or all of datasets. We provide Random Forests26, XGBoost27 and LightGBM28 with its corresponding parameters for basic tuning. The data will be split into training and test (20%) datasets by default. AUC-ROC29 or MSE30 will be shown according to different cases. This module will provide SHAP31 feature importance plot and SHAP dependence plot.

##### E – Double Machine Learning (DML) Analysis
We used DML introduced from RNAcompare to bypass data hamalisation for comparsion between different batches/tissues/omics levels. Causal Forests32 is recommended because it can estimate conditional treatment heterogeneity between batches, which can allow us to explore whether the most contributing factor (hub metapathways) in a given condition is the most important compared to protective mechanisms. Notably, the most significant contributor identified by Association analysis may differ from feature importance rankings obtained using Causal Forests.

##### F - Data Visualisation II
In this section, user can use network plot to visualise the interactions among metagenes/cell types for finding the most frequent/significant metagenes/cell types.

## Result & Demonstration
### Data Collection
#### 1. First we need to import ICA processed 1st level data from RNAcompare or user can calculate manually. ID_REF is the compulsory field.
   ![image](https://github.com/user-attachments/assets/f6c7700b-76db-4550-a616-21404425c389)
#### 2. Next we import the clinical data.
![image](https://github.com/user-attachments/assets/5a3dc8a9-5361-449c-b4e8-6df9aa6f0f18)

#### 3. Activate the labels apart from LABEL itself
![image](https://github.com/user-attachments/assets/50a9646e-763f-4f40-9bfa-e295e91a9cda)

#### 4. Run 2nd level ICA(metapathway)
![image](https://github.com/user-attachments/assets/b3adeda6-f804-4ffe-bb93-d8c840842abf)

and we can download the result.
![image](https://github.com/user-attachments/assets/56e95fc9-b6bb-4b20-b357-04d0d6461fc1)

#### 5. Data Visualisation I
We then can do the metapathway analysis to check the components of each and adjust the threshold.
![image](https://github.com/user-attachments/assets/01859e63-7efb-4884-87ee-2dc1b87efc10)

![image](https://github.com/user-attachments/assets/b8b57f70-f22d-4e84-998e-8e8f0fecd47a)

Here is the heatmap plot checking the relationship of metapathways.
![image](https://github.com/user-attachments/assets/1016fd2a-4e71-46d9-a0ba-b424f5ffd5d3)

Then, we can draw the matrixplot by each metapathways or each batch. In this example we only have one batch. So we draw it by metapathway.
![image](https://github.com/user-attachments/assets/4058350a-f78c-4847-b9a4-2bcaf67e72cc)

![image](https://github.com/user-attachments/assets/b0d3312f-6a05-4404-88b1-28b4c5662821)













