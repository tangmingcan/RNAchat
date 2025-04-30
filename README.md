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

## Result & Demonstration (Case I -drug response)
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

![image](https://github.com/user-attachments/assets/9ba46a04-5d05-4f9a-adae-1871e0859d66)

**Note, we already know that metagene_3 is IFN signaling from RNAcompare.**

Here is the heatmap plot checking the relationship of metapathways.
![image](https://github.com/user-attachments/assets/1016fd2a-4e71-46d9-a0ba-b424f5ffd5d3)

Then, we can draw the matrixplot by each metapathways or each batch. In this example we only have one batch. So we draw it by metapathway.

![image](https://github.com/user-attachments/assets/4058350a-f78c-4847-b9a4-2bcaf67e72cc)

![image](https://github.com/user-attachments/assets/b0d3312f-6a05-4404-88b1-28b4c5662821)

#### 6. Run Association analysis and check where is metapathway_5 (containing metagene_3) using SHAP plot.

![image](https://github.com/user-attachments/assets/ca80caec-1e6f-4cca-8ef2-e872f2f4282f)


![image](https://github.com/user-attachments/assets/9d76eff6-164a-4ab8-8bad-87af558bbb70)

#### 7. Now, let's draw Rituximab and anti-TNF separately. We use DML tab and treat them as 2 batches separately.

![image](https://github.com/user-attachments/assets/91c5420f-01e4-46fc-ba09-db119e021ba8)

Now, let's draw Rituximab.

![image](https://github.com/user-attachments/assets/94ca53a5-f6f6-48fb-9f51-72c0ab8087a1)

As we can see metapathway_2 is ordered before metapathway_5 now.

Let's draw Anti-TNF.

![image](https://github.com/user-attachments/assets/97b6fd4a-fd68-4ae6-909b-97b7b22ee0a3)

Here, matapathway_5 is ordered before metapathway_2.


The following is enrichment analysis  of metagene_3, metagene_2 from RNAcompare

**metagene_3**
![image](https://github.com/user-attachments/assets/35fd6a73-bbbf-4ad1-bf75-3b97bf67da8e)

**metagene_2**
![image](https://github.com/user-attachments/assets/e2c93dff-fa2a-4d88-b20d-1599df567f3b)

#### 8. Data Visualisation II- Interaction of metagenes
![image](https://github.com/user-attachments/assets/5204a5e2-f523-4988-901c-56e17173117c)


## Result & Demonstration (Case II -host vs parasite)
### Data collection
#### 1. We first upload 1st level ICA processed data for both plasmodium & human
![image](https://github.com/user-attachments/assets/0f258d2d-4cf1-4747-9779-3f7aafa0ad3a)

#### 2. We upload clinic data, and set 0 to the columns of parasites but keep the IDs.
![image](https://github.com/user-attachments/assets/99a1e862-be9b-4361-9a94-a7d7bca9e4e5)

#### 3. We ran metapathway, the 2nd level ICA.
![image](https://github.com/user-attachments/assets/0ef568a2-1c79-4e23-9c9a-72df8f8d17a1)

#### 4. We check the components of metapathways.
![image](https://github.com/user-attachments/assets/88ad5edd-0c9d-4603-aeed-295d33a90f72)

#### 5. We run association analysis using DML, because we only care about human batch.
![image](https://github.com/user-attachments/assets/2ae46b0f-33f5-4900-957e-3a393e51a103)

#### 6. We get SHAP plot.
![image](https://github.com/user-attachments/assets/c2afc19a-5efb-427f-ac88-ffbda8f845ad)

#### 7. Relationship plot.
![image](https://github.com/user-attachments/assets/a71817c2-b53a-4c6b-977b-4d5a0792f937)

![image](https://github.com/user-attachments/assets/88f9d341-d23b-46f8-a724-7e6a849b0e2c)

## Result & Demonstration (Case III - RA & MS)
We use a new case based on the paper [Unveiling the ageing-related genes in diagnosing osteoarthritis with metabolic syndrome by integrated bioinformatics analysis and machine learning](https://www.tandfonline.com/doi/full/10.1080/21691401.2025.2471762?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org#abstract), what we do will not only find the related signatures mentioned by the paper, but also we can find their corresponding pathways and metapathways in a supervised manner (importance to age).

The paper still has limitations: (1) It didn't fully use the clinical feature 'age' in the datasets, which can be used in the supervised ML rather than using an external aging-related Gene dataset to find overlapping genes. (2) Didn't find the corresponding pathways and metapathways and didn't prioritize them (not full transparent in terms of mechanism interpretation).(3) It needs control/healthy records for DE. RNAchat doesn't need healthy records!

#### 1. Data Upload for RNAcompare.
A slightly different from the datasets used by the paper, we used PEAC(synovium, refed by RNAcompare) and GSE58795(Metabolistic syndrom-PBMC)
**Note, PEAC is normalized(RNA-seq) and GSE58795 is raw data(microarray), but since we don't need batch correction any more, we can combine them via RNAcompare directly!(exciting?!)**

These are the data we load into RNAcompare:
![image](https://github.com/user-attachments/assets/6f6cd858-cbb0-4adf-8dea-eb4c98e57979)

![image](https://github.com/user-attachments/assets/6b496a60-bccc-43f3-9c19-22d15d962d7d)

#### 2. Data Upload for RNAchat.
We then do the 1st level ICA in RNAcompare, extracting 10 components for each disease. Then import this data to RNAchat.

![image](https://github.com/user-attachments/assets/87c4782d-3807-4bc8-a62d-05e764c8b5f6)

![image](https://github.com/user-attachments/assets/587ff89b-cfe9-4a94-9d3a-82e346cb2662)

#### 3. We ran metapathway, the 2nd level ICA and check the components.
![image](https://github.com/user-attachments/assets/4039cee5-f95e-4058-a298-00e6248285c6)

#### 4. We ran Causal Forests in terms of age.
![image](https://github.com/user-attachments/assets/3cd85f7b-607d-4eaa-851c-28ae662ca56c)

#### 5. We checked the SHAP importance plot.
Here we just care about the width of each feature. Then we can prioritize the importance of metagenes.
![image](https://github.com/user-attachments/assets/e082d083-ddfc-4d40-9fa3-1b01fb3cb8a2)

#### 6. Enrichment analysis for metagenes.
Metagene_8_ms:

![image](https://github.com/user-attachments/assets/b1231cae-ea33-4941-a74e-14caf9916d93)

You can check the important signatures whether they are in the supp material 1 of the original paper.

For the most important gene the original paper mentioned: CEBPB, we found it via the visualization plot:

![image](https://github.com/user-attachments/assets/be74280f-123a-4277-ab86-951848880040)

Here we can see metagene_4_RA is the most freq metagene interacting with others. So we checked its enrichment result.
You can see CEBPB is involved in pathway: WP_NETWORK_MAP_OF_SARSCOV2_SIGNALING
![image](https://github.com/user-attachments/assets/3dd8bb54-51fe-4ae2-bae1-ff8069785f26)














)

























